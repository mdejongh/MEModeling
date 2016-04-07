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
    use Data::Dumper;
    print &Dumper($result);
    done_testing(0);
};
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
