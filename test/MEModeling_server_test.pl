use strict;
use Data::Dumper;
use Test::More;
use Config::Simple;
use Time::HiRes qw(time);
use Bio::KBase::AuthToken;
use Bio::KBase::workspace::Client;
use MEModeling::MEModelingImpl;
use JSON::XS;

local $| = 1;
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('MEModeling');
my $ws_url = $config->{"workspace-url"};
my $ws_name = undef;
my $ws_client = new Bio::KBase::workspace::Client($ws_url,token => $token);
my $auth_token = Bio::KBase::AuthToken->new(token => $token, ignore_authrc => 1);
my $ctx = LocalCallContext->new($token, $auth_token->user_id);
$MEModeling::MEModelingServer::CallContext = $ctx;
my $impl = new MEModeling::MEModelingImpl();

sub get_ws_name {
    if (!defined($ws_name)) {
        my $suffix = int(time * 1000);
        $ws_name = 'test_MEModeling_' . $suffix;
        $ws_client->create_workspace({workspace => $ws_name});
    }
    return $ws_name;
}

eval {
    $@ = '';
    my $result;
    eval { 
	print STDERR "Loading genome and contigs ...\n";
	
        my $obj_name = "kb|g.0.c.1";
	open (CONTIG, "kb_g.0.contigset.json");
        my $obj = <CONTIG>;
	chomp $obj;
	close CONTIG;
	my $decoded = JSON::XS::decode_json($obj);
        $ws_client->save_objects({'workspace' => get_ws_name(), 'objects' => [{'type' => 'KBaseGenomes.ContigSet', 'name' => $obj_name, 'data' => $decoded}]});

	my $obj_name2 = "E_coli_K12_reannotated";
	open (ECK12, "E_coli_K12_reannotated.json");
        my $obj2 = <ECK12>;
	chomp $obj2;
	close ECK12;
	my $decoded2 = JSON::XS::decode_json($obj2);
	$decoded2->{"contigset_ref"} = get_ws_name()."/".$obj_name;
        $ws_client->save_objects({'workspace' => get_ws_name(), 'objects' => [{'type' => 'KBaseGenomes.Genome', 'name' => $obj_name2, 'data' => $decoded2}]});

	print STDERR "Getting ready ...\n";
        $result = $impl->build_me_model(get_ws_name, $obj_name2);
	print STDERR "Done\n";
    };
    print "$@\n";

    # apply petri-net test
    my %b0014 = petri('kb|g.0.peg.839',$result->{reactions},$result->{compounds});
    my @zeros = grep { $b0014{$_} == 0 } keys %b0014;
    delete @b0014{@zeros};
    print STDERR &Dumper(\%b0014);

    done_testing(0);
};

sub petri {
    my %free = map { $_ => 1 } qw/gtp 
Ala
Arg
Asn
Asp
Cys
Gln
Glu
Gly
His
Ile
Leu
Lys
Met
Phe
Pro
SeC(p)
Ser
Thr
Trp
Tyr
Val
 /;

    my ($gene,$reactions,$componds) = @_;
    $gene =~ s/\W/_/g;

    my %rxns2substrates;
    my %rxns2products;
    my %revrxns;
    my %recycling;

    sub myrib {
	if ($a =~ /tl_ini_.*_(\d+)_rib/) {
	    my $arib = $1;
	    if ($b =~ /tl_ini_.*_(\d+)_rib/) {
		my $brib = $1;
		return $arib <=> $brib;
	    }
	}

	return $a cmp $b;
    }

    my @reactions = (@{$reactions->{$gene}}, @{$reactions->{Recycling}});

    foreach (@reactions) {
	chomp;
	my ($id, $def, $equation, $rev, $subsys) = split "\t";
	my ($substrates, $products) = split /-->/, $equation;

	if ($subsys eq "Sinks" && $products eq "") {
	    $products = $substrates;
	    $substrates = "1 START_TOKEN";
	    $rev = "irreversible";
	}
	elsif ($subsys =~ "Recycling") {
	    $recycling{$id."_f"} = 1;
	    $recycling{$id."_r"} = 1;
	}

	while ($substrates =~ /(\d+)\s(\S+)/g) {
	    $rxns2substrates{$id."_f"}{$2} = $1;
	    if ($rev eq "reversible") {
		$rxns2products{$id."_r"}{$2} = $1;	    
		$revrxns{$id."_r"} = $id."_f";
		$revrxns{$id."_f"} = $id."_r";
	    }
	    elsif ($rev ne "irreversible") {
		print STDERR "strange: $rev\n";
	    }
	}
	while ($products =~ /(\d+)\s(\S+)/g) {
	    $rxns2products{$id."_f"}{$2} = $1;
	    if ($rev eq "reversible") {
		$rxns2substrates{$id."_r"}{$2} = $1;	    
	    }
	}
    }

    my %pool = ( "START_TOKEN" => 1, "EF-Ts" => 1 );
    my %used_rxns;
    my $recycle = 0;

    while (1) {
	my $fired = 0;

	foreach my $rxn (sort myrib keys %rxns2substrates) {
	    next if $recycle == 0 && $used_rxns{$rxn} == 1;
	    next if $recycle == 0 && exists $recycling{$rxn};
	    next if $recycle == 1 && ! exists $recycling{$rxn};
	    my $ready = 1;

	    foreach my $cpd (keys $rxns2substrates{$rxn}) {
		if ((($recycle == 0 && $cpd =~ /$gene/) || ($recycle == 1 && (! exists $free{$cpd}))) && $pool{$cpd} < $rxns2substrates{$rxn}{$cpd}) {
		    $ready = 0;
		    last;
		}
	    }

	    if ($ready == 1) {
		print STDERR "Firing $rxn\n" if $recycle == 0;
		foreach my $cpd (keys $rxns2substrates{$rxn}) {
		    $pool{$cpd} -= $rxns2substrates{$rxn}{$cpd};
		}
		foreach my $cpd (keys $rxns2products{$rxn}) {
		    $pool{$cpd} += $rxns2products{$rxn}{$cpd};
		}
		$fired = 1;
		$used_rxns{$rxn} = 1;
		$used_rxns{$revrxns{$rxn}} = 1;
	    }
	}

	if ($fired == 0) {
	    if ($recycle == 0) {
		print STDERR "Switching to recycling\n";
		$recycle = 1;
	    }
	    else {
		last;
	    }
	}
    }
    return %pool;
}

my $err = undef;
if ($@) {
    $err = $@;
}
eval {
    if (defined($ws_name)) {
        $ws_client->delete_workspace({workspace => $ws_name});
        print("Test workspace was deleted\n");
    }
};
if (defined($err)) {
    if(ref($err) eq "Bio::KBase::Exceptions::KBaseException") {
        die("Error while running tests: " . $err->trace->as_string);
    } else {
        die $err;
    }
}

{
    package LocalCallContext;
    use strict;
    sub new {
        my($class,$token,$user) = @_;
        my $self = {
            token => $token,
            user_id => $user
        };
        return bless $self, $class;
    }
    sub user_id {
        my($self) = @_;
        return $self->{user_id};
    }
    sub token {
        my($self) = @_;
        return $self->{token};
    }
    sub provenance {
        my($self) = @_;
        return [{'service' => 'MEModeling', 'method' => 'please_never_use_it_in_production', 'method_params' => []}];
    }
    sub authenticated {
        return 1;
    }
    sub log_debug {
        my($self,$msg) = @_;
        print STDERR $msg."\n";
    }
    sub log_info {
        my($self,$msg) = @_;
        print STDERR $msg."\n";
    }
}
