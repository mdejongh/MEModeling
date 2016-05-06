package MEModeling::MEModelingImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

MEModeling

=head1 DESCRIPTION

A KBase module: MEModeling

=cut

#BEGIN_HEADER
use Bio::KBase::AuthToken;
use Bio::KBase::workspace::Client;
use Config::IniFiles;
use Data::Dumper;

sub get_dna {
    my ($feature, $contigset) = @_;
    my $ret = $feature->{dna_sequence};
    if (! defined $ret || length $ret == 0) {
	my @ret;
	foreach my $loc (@{$feature->{location}}) {
	    my ($contig_id, $start, $strand, $len) = @$loc;
	    my $end;
	    if ($strand eq "+") {
		foreach my $contig (@{$contigset->{contigs}}) {
		    if ($contig->{id} eq $contig_id) {
			push @ret, substr($contig->{sequence}, ($start-1), $len);
		    }
		}
	    }
	    elsif ($strand eq "-") {
		foreach my $contig (@{$contigset->{contigs}}) {
		    if ($contig->{id} eq $contig_id) {
			my $tmp = scalar reverse substr($contig->{sequence}, ($start-$len), $len);
			$tmp =~ tr/ACGTacgt/TGCAtgca/;
			push  @ret, $tmp;
		    }
		}
	    }
	}
	$ret = join "", @ret;
    }
    return uc $ret;
}

sub avg_rna_count {
    my @cds = @_;
    my %rna_count;
    
    foreach my $cds (@cds) {
	my $cdsA = $cds =~ tr/A//;
	my $cdsC = $cds =~ tr/C//;
	my $cdsG = $cds =~ tr/G//;
	my $cdsU = $cds =~ tr/T//;
	my $sum = length $cds;

	$rna_count{'C'}+=($cdsG*10 + $cdsC*9 + $cdsA*10 + $cdsU*10);	
	$rna_count{'H'}+=($cdsG*11 + $cdsC*11 + $cdsA*11 + $cdsU*12);
	$rna_count{'N'}+=($cdsG*5 + $cdsC*3 + $cdsA*5 + $cdsU*2);
	$rna_count{'O'}+=($cdsG*7 + $cdsC*7 + $cdsA*6 + $cdsU*8-$sum);
	$rna_count{'P'}+=$sum;
	$rna_count{'charge'}+=$sum*(-1);
    }

    $rna_count{'C'} = int($rna_count{'C'}/@cds);
    $rna_count{'H'} = int($rna_count{'H'}/@cds);
    $rna_count{'N'} = int($rna_count{'N'}/@cds);
    $rna_count{'O'} = int($rna_count{'O'}/@cds);
    $rna_count{'P'} = int($rna_count{'P'}/@cds);
    $rna_count{'charge'} = int($rna_count{'charge'}/@cds);

    return \%rna_count;
}

sub parse_formula {
    my $formula = shift;
    my %retval;
    $retval{C}=0;
    $retval{H}=0;
    $retval{N}=0;
    $retval{O}=0;
    $retval{S}=0;
    $retval{P}=0;
    $retval{Mg}=0;
    $retval{Zn}=0;
    $retval{Fe}=0;
    my @f=split//,$formula;
    for (my $i=0;$i<=$#f;++$i)
    {
	if ($f[$i] =~ /\D/) 
	{
	    if ($i ne $#f)
	    {
		if ($f[$i+1] =~ /\D/)
		{
		    my $tmp=$f[$i];
		    $tmp.=$f[$i+1];
		    my $a='';
		    for (my $j=2;$j<=($#f-$i);++$j)
		    {
			if ($f[$i+$j] =~ /\d/)
			{
			    $a.=$f[$i+$j];
			}
			else
			{
			    $j=$#f;
			}
		    }
		    $retval{$tmp}=$a;
		    ++$i;
		    
		}
		else
		{	
		    $a='';
		    for (my $j=1;$j<=($#f-$i);++$j)
		    {
			if ($f[$i+$j] =~ /\d/)
			{
			    $a.=$f[$i+$j];
			}
			else
			{
			    $j=$#f;
			}
		    }
		    $retval{$f[$i]}=$a;
		}
	    }
	    else
	    {
		$retval{$f[$i]}=1;
	    }
	}
    }
    return \%retval;
}

sub aa_count {
    my ($aa_seq, $infos) = @_;
    my %infos = %$infos;
    my %aa_count;
    $aa_count{C}=0;
    $aa_count{H}=0;
    $aa_count{N}=0;
    $aa_count{O}=0;
    $aa_count{S}=0;
    $aa_count{Se}=0;
    $aa_count{charge}=0;

    foreach my $aas (keys %infos)
    {
	my $count = () = $aa_seq =~ /$aas/g;
	$aa_count{C}=$aa_count{C} + $infos{$aas}{C}*$count;
	$aa_count{H}=$aa_count{H} + $infos{$aas}{H}*$count;
	$aa_count{N}=$aa_count{N} + $infos{$aas}{N}*$count;
	$aa_count{O}=$aa_count{O} + $infos{$aas}{O}*$count;
	if (exists $infos{$aas}{S})
	{
	    $aa_count{S}=$aa_count{S} + $infos{$aas}{S}*$count;
	}
	if (exists $infos{$aas}{Se} && $infos{$aas}{Se} ne "")
	{
	    $aa_count{Se}=$aa_count{Se} + $infos{$aas}{Se}*$count;
	}
	$aa_count{charge}=$aa_count{charge} + $infos{$aas}{charge}*$count;
    }

    my $sum_aa = (length $aa_seq)-1; # subtracting one because ignoring first codon (fmet_tRNA or something like that)

    #adds formyl to AA formulae
    ++$aa_count{C};
    $aa_count{charge}=$aa_count{charge}-1;
    $aa_count{O}=($aa_count{O}-$sum_aa)+1+1;
    $aa_count{H}=$aa_count{H}-2*($sum_aa-1)-1;

    # comment from06-01-24: remove 1 h2o that is still present in formulae since $sum_aa does not account for 1st met
    $aa_count{H}=$aa_count{H}-2;
    $aa_count{O}=$aa_count{O}-1;

    return %aa_count;
}

sub mature_count {
    my %aa_count = @_;
    my %m_count;
    $m_count{C}=$aa_count{C}-1;
    $m_count{H}=$aa_count{H}+1;
    $m_count{N}=$aa_count{N};
    $m_count{O}=$aa_count{O}-1;
    $m_count{S}=$aa_count{S};
    $m_count{Se}=$aa_count{Se} if exists $aa_count{Se};
    $m_count{P}=0;
    $m_count{Mg}=0;
    $m_count{Fe}=0;
    $m_count{Zn}=0;
    $m_count{charge}=$aa_count{charge}+1;
    return %m_count;
}

#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR
    
    my $config_file = $ENV{ KB_DEPLOYMENT_CONFIG };
    my $cfg = Config::IniFiles->new(-file=>$config_file);
    my $wsInstance = $cfg->val('MEModeling','workspace-url');
    die "no workspace-url defined" unless $wsInstance;
    
    $self->{'workspace-url'} = $wsInstance;
    
    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 build_me_model

  $return = $obj->build_me_model($workspace, $genome_id)

=over 4

=item Parameter and return types

=begin html

<pre>
$workspace is a MEModeling.workspace_name
$genome_id is a MEModeling.ws_genome_id
$return is a MEModeling.MEModelingResult
workspace_name is a string
ws_genome_id is a string
MEModelingResult is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string

</pre>

=end html

=begin text

$workspace is a MEModeling.workspace_name
$genome_id is a MEModeling.ws_genome_id
$return is a MEModeling.MEModelingResult
workspace_name is a string
ws_genome_id is a string
MEModelingResult is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string


=end text



=item Description

Build a Metabolic/Expression model

=back

=cut

sub build_me_model
{
    my $self = shift;
    my($workspace, $genome_id) = @_;

    my @_bad_arguments;
    (!ref($workspace)) or push(@_bad_arguments, "Invalid type for argument \"workspace\" (value was \"$workspace\")");
    (!ref($genome_id)) or push(@_bad_arguments, "Invalid type for argument \"genome_id\" (value was \"$genome_id\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to build_me_model:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'build_me_model');
    }

    my $ctx = $MEModeling::MEModelingServer::CallContext;
    my($return);
    #BEGIN build_me_model

    my $token=$ctx->token;
    my $wsClient=Bio::KBase::workspace::Client->new($self->{'workspace-url'},token=>$token);
    my $provenance = [{}];
    $provenance = $ctx->provenance if defined $ctx->provenance;
    
    my $genome;
    eval {
	$genome = $wsClient->get_objects([{ref=>$workspace."/".$genome_id}])->[0]{data};
	push @{$provenance->[0]->{'input_ws_objects'}}, $workspace."/".$genome_id;
    };
    if ($@) {
	die "Error loading genome:\n".$@;
    }

    my $contigset;
    eval {
	$contigset = $wsClient->get_objects([{ref=>$genome->{contigset_ref}}])->[0]{data};
    };
    if ($@) {
	die "Error loading contigset:\n".$@;
    }

    # map compound abbrevs to modelseed ids
    my %cpd_map = (
	  'mg2' => 'cpd00254',
          'amp' => 'cpd00018',
          'Lys' => 'cpd00039',
          'Arg' => 'cpd00051',
          'Glu' => 'cpd00023',
          'Asp' => 'cpd00041',
          'Phe' => 'cpd00066',
          'Met' => 'cpd00060',
          'Ile' => 'cpd00322',
          'Val' => 'cpd00156',
          'Gly' => 'cpd00033',
          'His' => 'cpd00119',
          'gtp' => 'cpd00038',
          'gdp' => 'cpd00031',
          'ppi' => 'cpd00012',
          'ump' => 'cpd00091',
          'Tyr' => 'cpd00069',
          'Asn' => 'cpd00132',
          'Pro' => 'cpd00129',
          'atp' => 'cpd00002',
          'gmp' => 'cpd00126',
          'utp' => 'cpd00062',
          'Thr' => 'cpd00161',
          'Ala' => 'cpd00035',
          'pi' => 'cpd00009',
          'Ser' => 'cpd00054',
          'h2o' => 'cpd00001',
          'adp' => 'cpd00008',
          'Gln' => 'cpd00053',
          'ctp' => 'cpd00052',
          'h' => 'cpd00067',
          'Leu' => 'cpd00107',
          'cmp' => 'cpd00046',
          'Cys' => 'cpd00084',
          'Trp' => 'cpd00065'
        );

    # NEED GENETIC CODE
    my $assigned_code = $genome->{genetic_code};

    ## Need FR to GENE so we can calculate formulae for factors
    ## Also need dna sequences for tRNAs and protein sequences for pegs
    my %aa_seq;
    my %fr2gene;
    my %trna_seqs;

    foreach my $feature (@{$genome->{features}}) {
	my $gene = $feature->{id};
	$gene =~ s/\W/_/g;
	my $fr = $feature->{function};
	my $dna = &get_dna($feature, $contigset);
	$fr =~ s/\s*#.*$//;
	if ($fr =~ /\s*;\s*/) {
	    map { $fr2gene{$_} = $gene } split /s*;\s*/, $fr;
	}
	elsif ($fr =~ /\s*@\s*/) {
	    map { $fr2gene{$_} = $gene } split /\s*@\s*/, $fr;
	}
	else {
	    $fr2gene{$fr} = $gene;
	}
	if ($fr =~ /^tRNA\-(\w+)\-\w{3}$/) {
	    push @{$trna_seqs{$fr}}, uc $dna;
	}
	else {
	    $aa_seq{$gene} = uc $feature->{protein_translation};
	}

    }

    # Compute average composition of tRNAs that match the same codons 
    my %formula_tRNA;

    foreach my $fr (keys %trna_seqs) {
	$formula_tRNA{$fr} = &avg_rna_count(@{$trna_seqs{$fr}});
    }

    ###########################################
    # Amino Acid Info
    ###########################################
    my @AAInfo;
    $AAInfo[1] = 'A,ala-L,C3H7N1O2,0';
    $AAInfo[2] = 'R,arg-L,C6H15N4O2,1';
    $AAInfo[3] = 'N,asn-L,C4H8N2O3,0';
    $AAInfo[4] = 'D,asp-L,C4H6N1O4,-1';
    $AAInfo[5] = 'C,cys-L,C3H7N1O2S1,0';
    $AAInfo[6] = 'Q,gln-L,C5H10N2O3,0';
    $AAInfo[7] = 'E,glu-L,C5H8N1O4,-1';
    $AAInfo[8] = 'G,gly,C2H5N1O2,0';
    $AAInfo[9] = 'H,his-L,C6H9N3O2,0';
    $AAInfo[10] = 'I,ile-L,C6H13N1O2,0';
    $AAInfo[11] = 'L,leu-L,C6H13N1O2,0';
    $AAInfo[12] = 'K,lys-L,C6H15N2O2,1';
    $AAInfo[13] = 'M,met-L,C5H11N1O2S1,0';
    $AAInfo[14] = 'F,phe-L,C9H11N1O2,0';
    $AAInfo[15] = 'P,pro-L,C5H9N1O2,0';
    $AAInfo[16] = 'S,ser-L,C3H7N1O3,0';
    $AAInfo[17] = 'T,thr-L,C4H9N1O3,0';
    $AAInfo[18] = 'W,trp-L,C11H12N2O2,0';
    $AAInfo[19] = 'Y,tyr-L,C9H11N1O3,0';
    $AAInfo[20] = 'V,val-L,C5H11N1O2,0';
    $AAInfo[21] = 'Sec,Sec,C3H7N1O2Se1,0';

    my %infos;
    for (my $k=1;$k<=$#AAInfo;++$k)
    {
	my @data=split ",", $AAInfo[$k];
	$infos{$data[0]}=&parse_formula($data[2]);
	$infos{$data[0]}{charge}=$data[3];
    }

    my %aa_abbrev;
    $aa_abbrev{"Ala"} = "A";
    $aa_abbrev{"Arg"} = "R";
    $aa_abbrev{"Asn"} = "N";
    $aa_abbrev{"Asp"} = "D";
    $aa_abbrev{"Cys"} = "C";
    $aa_abbrev{"Gln"} = "Q";
    $aa_abbrev{"Glu"} = "E";
    $aa_abbrev{"Gly"} = "G";
    $aa_abbrev{"His"} = "H";
    $aa_abbrev{"Ile"} = "I";
    $aa_abbrev{"Leu"} = "L";
    $aa_abbrev{"Lys"} = "K";
    $aa_abbrev{"Met"} = "M";
    $aa_abbrev{"Phe"} = "F";
    $aa_abbrev{"Pro"} = "P";
    $aa_abbrev{"SeC(p)"} = "Sec";
    $aa_abbrev{"Ser"} = "S";
    $aa_abbrev{"Thr"} = "T";
    $aa_abbrev{"Trp"} = "W";
    $aa_abbrev{"Tyr"} = "Y";
    $aa_abbrev{"Val"} = "V";

    # read factors.txt and calculate formula and charge
    my (%factors, %formulae);

    open (FACTORS, "../data/factors.txt") or die("Couldn't open factors.txt: $!");
    while (<FACTORS>) {
	chomp;
	my ($category, $fname, $fr, $formula, $charge) = split "\t";
	if ($fr ne "" && ! exists $fr2gene{$fr}) {
	    print STDERR "No gene defined for TT factor '$fr' [$fname]\n";
	    next;
	}
	my ($extra, %extra);
	if ($fname =~ /\.GDP$/) {
	    $extra = "_GDP";
	    %extra = ( 'C' => 10, 'H' => 13, 'N' => 5, 'O' => 11, 'P' => 2, 'charge' => -2 );
	}
	elsif ($fname =~ /\.GTP$/) {
	    $extra = "_GTP";
	    %extra = ( 'C' => 10, 'H' => 16, 'N' => 5, 'O' => 14, 'P' => 3, 'charge' => -3 );
	}
	elsif ($fname =~ /\.me\.spmd$/) {
	    $extra = "_me_spmd";
	    %extra = ( 'C' => 11, 'H' => 33, 'N' => 4, 'charge' => 3 );
	}
	elsif ($fname =~ /\.spmd$/) {
	    $extra = "_spmd";
	    %extra = ( 'C' => 10, 'H' => 30, 'N' => 4, 'charge' => 4 );
	}
	if (defined $category) {
	    if ($fr ne "") {
		$factors{$category}{$fname} = $fr2gene{$fr}.$extra;
	    }
	    else {
		$factors{$category}{$fname} = $fname;
	    }
	}
	if (! defined $formula) {
	    my %m_count = &mature_count(&aa_count($aa_seq{$fr2gene{$fr}}, \%infos));
	    my $multiplier = 1;
	    if ($fname =~ /_hexa$/) {
		$multiplier = 6;
	    }
	    elsif ($fname =~ /_dim$/) {
		$multiplier = 2;
	    }
	    if ($multiplier > 1) {
		foreach my $key (keys %m_count) {
		    $m_count{$key} *= $multiplier;
		}
	    }
	    if (keys %extra > 0) {
		foreach my $key (keys %extra) {
		    $m_count{$key} += $extra{$key};
		}
	    }
	    $formulae{$fname} = \%m_count;
	}
	else {
	    $formulae{$fname} = parse_formula($formula);
	    $formulae{$fname}{charge} = $charge;
	}
    }

    # now we can calculate the EF-Tu.GTP.trnas
    foreach my $trna (keys %formula_tRNA) {
	my %formula = %{$formulae{"EF-Tu.GTP"}};
	if ($trna =~ /^tRNA\-(\w+)\-\w{3}$/) {
	    my $abbrev = $aa_abbrev{$1};
	    map { $formula{$_} += $formula_tRNA{$trna}{$_} } keys %{$formula_tRNA{$trna}};
	    map { $formula{$_} += $infos{$abbrev}{$_} } keys %{$infos{$abbrev}};
	}
	$formulae{"EF-Tu.GTP.$trna"} = \%formula;
    }

    ###########################################
    # Molecular Weight of different Elements
    ###########################################
    my %mw;
    $mw{C}=12.0107;
    $mw{H}=1.00794;
    $mw{O}=15.9994;
    $mw{N}=14.0067;
    $mw{P}=30.973762;
    $mw{Mg}=24.3050;
    $mw{Zn}=65.409;
    $mw{Fe}=55.845;
    $mw{S}=32.065;
    $mw{Se}=78.96;

    ## GENETIC CODE
    my %genetic_code;
    map { $genetic_code{'11'}{'Ala'}{$_} = 1 } ('GCT','GCA','GCG','GCC');
    map { $genetic_code{'11'}{'Arg'}{$_} = 1 } ('CGT','CGC','CGA','AGA','AGG','CGG');
    map { $genetic_code{'11'}{'Asn'}{$_} = 1 } ('AAC','AAT');
    map { $genetic_code{'11'}{'Asp'}{$_} = 1 } ('GAC','GAT');
    map { $genetic_code{'11'}{'Cys'}{$_} = 1 } ('TGC','TGT');
    map { $genetic_code{'11'}{'Gln'}{$_} = 1 } ('CAG','CAA');
    map { $genetic_code{'11'}{'Glu'}{$_} = 1 } ('GAA','GAG');
    map { $genetic_code{'11'}{'Gly'}{$_} = 1 } ('GGC','GGT','GGA','GGG');
    map { $genetic_code{'11'}{'His'}{$_} = 1 } ('CAC','CAT');
    map { $genetic_code{'11'}{'Ile'}{$_} = 1 } ('ATC','ATT','ATA');
    map { $genetic_code{'11'}{'Leu'}{$_} = 1 } ('CTG','CTC','CTT','CTA','TTG','TTA');
    map { $genetic_code{'11'}{'Lys'}{$_} = 1 } ('AAA','AAG');
    map { $genetic_code{'11'}{'Met'}{$_} = 1 } ('ATG');
    map { $genetic_code{'11'}{'Phe'}{$_} = 1 } ('TTC','TTT');
    map { $genetic_code{'11'}{'Pro'}{$_} = 1 } ('CCG','CCC','CCT','CCA');
    map { $genetic_code{'11'}{'Ser'}{$_} = 1 } ('TCC','TCT','TCA','TCG','AGC','AGT');
    map { $genetic_code{'11'}{'Thr'}{$_} = 1 } ('ACC','ACT','ACA','ACG');
    map { $genetic_code{'11'}{'Trp'}{$_} = 1 } ('TGG');
    map { $genetic_code{'11'}{'Tyr'}{$_} = 1 } ('TAC','TAT');
    map { $genetic_code{'11'}{'Val'}{$_} = 1 } ('GTA','GTG','GTC','GTT');

    map { $genetic_code{'4'}{'Ala'}{$_} = 1 } ('GCT','GCA','GCG','GCC');
    map { $genetic_code{'4'}{'Arg'}{$_} = 1 } ('CGT','CGC','CGA','AGA','AGG','CGG');
    map { $genetic_code{'4'}{'Asn'}{$_} = 1 } ('AAC','AAT');
    map { $genetic_code{'4'}{'Asp'}{$_} = 1 } ('GAC','GAT');
    map { $genetic_code{'4'}{'Cys'}{$_} = 1 } ('TGC','TGT');
    map { $genetic_code{'4'}{'Gln'}{$_} = 1 } ('CAG','CAA');
    map { $genetic_code{'4'}{'Glu'}{$_} = 1 } ('GAA','GAG');
    map { $genetic_code{'4'}{'Gly'}{$_} = 1 } ('GGC','GGT','GGA','GGG');
    map { $genetic_code{'4'}{'His'}{$_} = 1 } ('CAC','CAT');
    map { $genetic_code{'4'}{'Ile'}{$_} = 1 } ('ATC','ATT','ATA');
    map { $genetic_code{'4'}{'Leu'}{$_} = 1 } ('CTG','CTC','CTT','CTA','TTG','TTA');
    map { $genetic_code{'4'}{'Lys'}{$_} = 1 } ('AAA','AAG');
    map { $genetic_code{'4'}{'Met'}{$_} = 1 } ('ATG');
    map { $genetic_code{'4'}{'Phe'}{$_} = 1 } ('TTC','TTT');
    map { $genetic_code{'4'}{'Pro'}{$_} = 1 } ('CCG','CCC','CCT','CCA');
    map { $genetic_code{'4'}{'Ser'}{$_} = 1 } ('TCC','TCT','TCA','TCG','AGC','AGT');
    map { $genetic_code{'4'}{'Thr'}{$_} = 1 } ('ACC','ACT','ACA','ACG');
    map { $genetic_code{'4'}{'Trp'}{$_} = 1 } ('TGG','TGA');
    map { $genetic_code{'4'}{'Tyr'}{$_} = 1 } ('TAC','TAT');
    map { $genetic_code{'4'}{'Val'}{$_} = 1 } ('GTA','GTG','GTC','GTT');

    my %LastCodon;
    $LastCodon{"TAA"}=1;
    $LastCodon{"TAG"}=1;
    my $rf1_mono = $factors{LastCodonsFactors}{"RF1_mono"};
    delete $factors{LastCodonsFactors}{"RF1_mono"};
    my $rf2_mono = $factors{LastCodonsFactors}{"RF2_mono"};
    delete $factors{LastCodonsFactors}{"RF2_mono"};
    $factors{LastCodonsFactors}{"RF1_mono"}{"TAA"}=$rf1_mono;
    $factors{LastCodonsFactors}{"RF1_mono"}{"TAG"}=$rf1_mono;
    $factors{LastCodonsFactors}{"RF2_mono"}{"TAA"}=$rf2_mono;

    if ($assigned_code == '11') {
	$factors{LastCodonsFactors}{"RF2_mono"}{"TGA"}=$rf2_mono;
	$LastCodon{"TGA"}=1;
    }

    ## BEGIN LOAD FUNCTIONAL ROLES for tRNAs
    my %list_tRNA;
    my %code_tRNA;

    # try Gary's rules
    my %cm = ( "A" => "T", "T" => "A", "G" => "C", "C" => "G" );
    my %default_codon_assignment; # just in case we leave a codon unassigned after applying Gary's rules

    foreach my $fr (keys %formula_tRNA) {
	if ($fr =~ /^tRNA-(.*)-([ACGT]{3})$/) {
	    my $aa_name = $1;
	    my $anti_codon = $2;
	    if (! exists $aa_abbrev{$aa_name}) {
		print STDERR "No abbreviation available for amino acid name: $aa_name\n";
	    }
	    else {
		$list_tRNA{$fr} = $aa_abbrev{$aa_name};

		my %gary_codons;
		if ($aa_name eq "Met" && $anti_codon eq "CAT") {
		    $gary_codons{"ATG"} = 1;
		}
		elsif ($aa_name eq "Ile" && $anti_codon eq "CAT") {
		    $gary_codons{"ATA"} = 1;
		}
		elsif ($aa_name eq "SeC(p)" && $anti_codon eq "TCA") {
		    $gary_codons{"AGT"} = 1;
		}
		elsif ($anti_codon =~ /^A(.)(.)/) {
		    $gary_codons{$cm{$2}.$cm{$1}."T"} = 1;
		}
		elsif ($anti_codon =~ /^C(.)(.)/) {
		    $gary_codons{$cm{$2}.$cm{$1}."G"} = 1;
		}
		elsif ($anti_codon =~ /^G(.)(.)/) {
		    $gary_codons{$cm{$2}.$cm{$1}."C"} = 1;
		    $gary_codons{$cm{$2}.$cm{$1}."T"} = 1;
		}
		elsif ($anti_codon =~ /^T(.)(.)/) {
		    $gary_codons{$cm{$2}.$cm{$1}."A"} = 1;
		    $gary_codons{$cm{$2}.$cm{$1}."G"} = 1;
		}

		# at this point we will only assign codons that are in the genetic code
		# AND are in gary's code
		foreach my $codon (keys %{$genetic_code{$assigned_code}{$aa_name}}) {
		    if (exists $gary_codons{$codon}) {
			$code_tRNA{$codon}{trna} = $fr;
			$code_tRNA{$codon}{aa} = $aa_abbrev{$aa_name};
		    }
		    else {
			# just in case it doesn't get assigned later
			$default_codon_assignment{$codon}{trna} = $fr;
			$default_codon_assignment{$codon}{aa} = $aa_abbrev{$aa_name};
		    }
		}
	    }
	}
    }

    # now make sure that there all codons are assigned
    foreach my $aa_name (keys %{$genetic_code{$assigned_code}}) {
	foreach my $codon (keys %{$genetic_code{$assigned_code}{$aa_name}}) {
	    if (! exists $code_tRNA{$codon}) {
		print STDERR "No tRNA assignment for codon $codon\n";
		if (exists $default_codon_assignment{$codon}) {
		    print STDERR "\t assigning $default_codon_assignment{$codon}{trna}\n";
		    $code_tRNA{$codon}{trna} = $default_codon_assignment{$codon}{trna};
		    $code_tRNA{$codon}{aa} = $default_codon_assignment{$codon}{aa};
		}
		else {
		    print STDERR "\t NO DEFAULT ASSIGNMENT!\n";
		}
	    }
	}
    }

    # sets sigma factor 70 as default sigma factor
    my $rnap='RNAP_70';
    my $sigm=$factors{TranscriptionSigmaFactor}{RpoD_mono};

    my (%reactions, %compounds);

    foreach my $feature (@{$genome->{features}}) {
	my $gene = $feature->{id};
	$gene =~ s/\W/_/g;
	my $fr = $feature->{function};
	my $type = $feature->{type};
	my $cds = &get_dna($feature, $contigset);
#	print STDERR "Processing $gene [$fr] of type $type\n";

	# CHECK THIS
	$type = "tRNA" if $type =~ /rna/ && $fr =~ /^tRNA-.*-[ACGT]{3}$/;

	# count occurrence of each type of base in CDS
	my $cdsA = $cds =~ tr/A//;
	my $cdsC = $cds =~ tr/C//;
	my $cdsG = $cds =~ tr/G//;
	my $cdsU = $cds =~ tr/T//;
	my $sum = length $cds;

	# count occurrence of each type of base in first 16 positions - IS A CHECK FOR DIRECTIONALITY NECESSARY?
	my $first16 = substr($cds,0,16);
	my $firstA = $first16 =~ tr/A//;
	my $firstC = $first16 =~ tr/C//;
	my $firstG = $first16 =~ tr/G//;
	my $firstU = $first16 =~ tr/T//;

	# DNA binding of activator
	my $cdscompC=$cdsG*10 + $cdsC*9 + $cdsA*10 + $cdsU*10;	
	my $cdscompH=$cdsG*11 + $cdsC*11 + $cdsA*11 + $cdsU*12; #i checked protons - 11-15
	my $cdscompN=$cdsG*5 + $cdsC*3 + $cdsA*5 + $cdsU*2;
	my $cdscompO=$cdsG*7 + $cdsC*7 + $cdsA*6 + $cdsU*8-$sum; # i changed this - 11-15
	my $cdscompP=$sum;
	my $cdscharge=$sum*(-1);
	push @{$compounds{$gene}}, "$gene\_DNA_act\tDNA transcription unit $gene (activated form)\tC${cdscompC}H${cdscompH}N${cdscompN}O${cdscompO}P${cdscompP}\t${cdscharge}\tTranscription\t, ${cdsG}G ${cdsC}C ${cdsU}U ${cdsA}A\n";
	push @{$compounds{$gene}}, "$gene\_DNA_neu\tDNA transcription unit $gene (inactive form)\tC${cdscompC}H${cdscompH}N${cdscompN}O${cdscompO}P${cdscompP}\t${cdscharge}\tTranscription\t, ${cdsG}G ${cdsC}C ${cdsU}U ${cdsA}A\n";
	push @{$reactions{$gene}}, "$gene\_DNA_act_bind\t$gene\_DNA binding of activator\t1 $gene\_DNA_neu --> 1 $gene\_DNA_act\tirreversible\tTranscription Regulation\n";
	push @{$reactions{$gene}}, "sink_$gene\_DNA_neu\tsink $gene\_DNA_neu\t1 $gene\_DNA_neu -->\treversible\tSinks\n";

	# initiation
	# complexes include frna and RNAP
	my $complexC = ${cdscompC}+${firstG}*10+ ${firstC}*9 + ${firstA} *10 + ${firstU} *9;
	my $complexH = ${cdscompH}+${firstG}*11 + ${firstC}*11 + ${firstA}  *11 + ${firstU} *10+1; 
	my $complexN = ${cdscompN}+${firstG}*5 + ${firstC}*3 + ${firstA}*5 + ${firstU}*2;
	my $complexO = ${cdscompO}+${firstG}*8 + ${firstC} * 8 + ${firstA}*7 + ${firstU}*9+6-15;
	my $complexS = 0;
	my $complexP = ${cdscompP}+18;
	my $complexMg = 0;
	my $complexZn = 0;
	my $complexFe = 0;
	my $complexcharge = ${cdscharge}-19;

	foreach my $factor (sort keys %{$factors{RNAP}})
	{
	    ${complexC}+=$formulae{$factor}{C};
	    ${complexH}+=$formulae{$factor}{H};
	    ${complexN}+=$formulae{$factor}{N};
	    ${complexO}+=$formulae{$factor}{O};
	    ${complexS}+=$formulae{$factor}{S};
	    ${complexP}+=$formulae{$factor}{P};
	    ${complexMg}+=$formulae{$factor}{Mg};
	    ${complexZn}+=$formulae{$factor}{Zn};
	    ${complexFe}+=$formulae{$factor}{Fe};
	    ${complexcharge}+=$formulae{$factor}{charge};
	}
	
	push @{$compounds{$gene}},"transcr_ini_$gene\_cplx\ttranscription initiation complex $gene\tC${complexC}H${complexH}N${complexN}O${complexO}S${complexS}P${complexP}Mg${complexMg}Zn${complexZn}Fe${complexFe}\t${complexcharge}\tTranscription\t${firstG}G ${firstC}C ${firstU}U ${firstA}A\n";

	if ($type eq "CDS") {
	    push @{$reactions{$gene}}, "tscr_ini_$gene\tTranscription initiation of $gene\t1 $gene\_DNA_act + 1 $rnap + $firstA $cpd_map{atp} + $firstC $cpd_map{ctp} + $firstG $cpd_map{gtp} + $firstU $cpd_map{utp} --> 1 transcr_ini_$gene\_cplx + 1 $sigm + ".($firstA + $firstC + $firstG + $firstU-1)." $cpd_map{ppi}\treversible\tTranscription\n";

	    ###########################
	    # RHO DEPENDENT TERMINATION
	    ###########################

	    # formation complex for elongation - start with initiation complex
	    my $NameFactors = '';			
	    foreach my $factor (sort keys %{$factors{RhoDependentTranscriptionTerminationCDS}})
	    {
		${complexC}+=$formulae{$factor}{C};
		${complexH}+=$formulae{$factor}{H};
		${complexN}+=$formulae{$factor}{N};
		${complexO}+=$formulae{$factor}{O};
		${complexS}+=$formulae{$factor}{S};
		${complexP}+=$formulae{$factor}{P};
		${complexMg}+=$formulae{$factor}{Mg};
		${complexZn}+=$formulae{$factor}{Zn};
		${complexFe}+=$formulae{$factor}{Fe};
		${complexcharge}+=$formulae{$factor}{charge};
		$NameFactors.= $factor.", ";
	    }	
	    push @{$compounds{$gene}}, "transcr_elo_$gene\_cplx\ttranscription elongation complex $gene (CDS, RHO DEPENDENT TERMINATION, $NameFactors)\tC${complexC}H${complexH}N${complexN}O${complexO}S${complexS}P${complexP}Mg${complexMg}Zn${complexZn}Fe${complexFe}\t${complexcharge}\tTranscription\t${firstG}G ${firstC}C ${firstU}U ${firstA}A\n";
	    my $rxn = "tscr_elo_$gene\_ini_rho_dep\tFormation complex for elongation of $gene (RHO DEPENDENT TERMINATION)\t1 transcr_ini_$gene\_cplx ";
	    foreach my $factor (sort keys %{$factors{RhoDependentTranscriptionTerminationCDS}})
	    {
		$rxn .= "+ 1 $factors{RhoDependentTranscriptionTerminationCDS}{$factor} ";
	    }		
	    $rxn .= "--> 1 transcr_elo_$gene\_cplx\treversible\tTranscription\n";
	    push @{$reactions{$gene}}, $rxn;

	    # elongation and termination
	    my $rxn2 = "tscr_elo_term_$gene\_rho_dep\tTranscription elongation and RHO DEPENDENT termination of $gene\t1 transcr_elo_$gene\_cplx + 3 $cpd_map{h2o} + ".($cdsA-$firstA+3)." $cpd_map{atp} + ".($cdsC-$firstC)." $cpd_map{ctp} + ".($cdsG-$firstG)." $cpd_map{gtp} + ".($cdsU-$firstU)." $cpd_map{utp} --> 1 $gene\_DNA_neu + 1 $gene\_mRNA + ".(($cdsA-$firstA)+($cdsC-$firstC)+($cdsG-$firstG)+($cdsU-$firstU))." $cpd_map{ppi} ";
	    foreach my $factor (sort keys %{$factors{RhoDependentTranscriptionTerminationCDS}})
	    {
		$rxn2 .= "+ 1 $factors{RhoDependentTranscriptionTerminationCDS}{$factor} ";
	    }		
	    foreach my $factor (sort keys %{$factors{RNAP}})
	    {
		$rxn2 .= "+ 1 $factors{RNAP}{$factor} ";
	    }
	    $rxn2 .= "+ 3 $cpd_map{adp} + 3 $cpd_map{pi} + 3 $cpd_map{h}\tirreversible\tTranscription\t".(($cdsA+$cdsC+$cdsG+$cdsU)/45)." (45nt/s) to ".(($cdsA+$cdsC+$cdsG+$cdsU)/40)." (40nt/s) 1/s\n";
	    push @{$reactions{$gene}}, $rxn2;
	}
	elsif ($type =~ /RNA/i)
	{

	    # assume this is true because we treat each gene as its own transcription unit: $gen_tmp{$cds_name}{first_gene } == 1
	    #***********************************
	    #RNA
	    #***********************************
	    my $sum_dna = ${cdsG} + ${cdsC}+${cdsU}+${cdsA};
	    my ${cds_rna_C}=${cdsG}*10 + ${cdsC}*9 + ${cdsA}*10 + ${cdsU}*9;	
	    my ${cds_rna_H}=${cdsG}*11 + ${cdsC}*11 + ${cdsA}*11 + ${cdsU}*10+1;
	    my ${cds_rna_N}=${cdsG}*5 + ${cdsC}*3 + ${cdsA}*5 + ${cdsU}*2;
	    my ${cds_rna_O}=${cdsG}*7 + ${cdsC}*7 + ${cdsA}*6 + ${cdsU}*8+7;
	    my ${cds_rna_P}=$sum_dna+2;
	    my ${cds_rna_charge}=$sum_dna*(-1)-3;
	    my $sum_rna=$sum_dna;
	    my $status_quo=''; # was 'first'
	    #${cds_rna_O}=${cds_rna_O}-6;
	    #${cds_rna_P}=${cds_rna_P}-2;
	    #${cds_rna_charge}=${cds_rna_charge}+2;

	    my $NameFactors = '';			
	    foreach my $factor (sort keys %{$factors{TranscriptionTerminationRNA}})
	    {
		${complexC}+=$formulae{$factor}{C};
		${complexH}+=$formulae{$factor}{H};
		${complexN}+=$formulae{$factor}{N};
		${complexO}+=$formulae{$factor}{O};
		${complexS}+=$formulae{$factor}{S};
		${complexP}+=$formulae{$factor}{P};
		${complexMg}+=$formulae{$factor}{Mg};
		${complexZn}+=$formulae{$factor}{Zn};
		${complexFe}+=$formulae{$factor}{Fe};
		${complexcharge}+=$formulae{$factor}{charge};
		$NameFactors.= $factor.", ";
	    }			
	    # transcr_elo_$gene_cplx accounts for: DNA + 16 mRNA, hRNAP, NusA_mono, NusG_mono, GreA_mono, GreB_mono, omega, Mfd_mono
	    push @{$compounds{$gene}},"transcr_elo_$gene\_cplx\ttranscription elongation complex ($NameFactors) $gene\tC${complexC}H${complexH}N${complexN}O${complexO}S${complexS}P${complexP}Mg${complexMg}Zn${complexZn}Fe${complexFe}\t${complexcharge}\tTranscription\t${firstG}G ${firstC}C ${firstU}U ${firstA}A\n";

	    push @{$compounds{$gene}}, "$gene\_RNA\t$type\tC${cds_rna_C}H${cds_rna_H}N${cds_rna_N}O${cds_rna_O}P${cds_rna_P}\t${cds_rna_charge}\tRNA cutting".(${cdsG}+${cdsC}+${cdsU}+${cdsA}).", ($status_quo), , ${cdsG}G ${cdsC}C ${cdsU}U ${cdsA}A\n";

	    push @{$reactions{$gene}}, "tscr_ini_$gene\_stab\tTranscription initiation of $gene (stable RNA) (, $gene)\t1 $gene\_DNA_act + 1 $rnap + ${firstA} $cpd_map{atp} + ${firstC} $cpd_map{ctp} + ${firstG} $cpd_map{gtp} + ${firstU} $cpd_map{utp} --> 1 transcr_ini_$gene\_cplx + 1 $sigm + ".(${firstA} + ${firstC}+ ${firstG} + ${firstU}-1)." $cpd_map{ppi}\treversible\tTranscription\n";
	    
	    my $rxn3 = "tscr_elo_$gene\_ini\_stab\tFormation complex for elongation of $gene (stable RNA)\t1 transcr_ini_$gene\_cplx ";
	    foreach my $factor (sort keys %{$factors{TranscriptionTerminationRNA}})
	    {
		$rxn3 .= "+ 1 $factors{TranscriptionTerminationRNA}{$factor} ";
	    }						
	    $rxn3 .= "--> 1 transcr_elo_$gene\_cplx\treversible\tTranscription\n";
	    push @{$reactions{$gene}}, $rxn3;
	    
	    my $rxn4 = "tscr_elo_term_$gene\_stab\tTranscription elongation and termination of $gene (stable RNA) (, $gene)\t1 transcr_elo_$gene\_cplx + ".(${cdsA}-${firstA})." $cpd_map{atp} + ".(${cdsC}-${firstC})." $cpd_map{ctp} + ".(${cdsG}-${firstG})." $cpd_map{gtp} + ".(${cdsU}-${firstU})." $cpd_map{utp} --> 1 $gene\_DNA_neu + 1 $gene\_RNA + ".((${cdsA}-${firstA})+(${cdsC}-${firstC})+(${cdsG}-${firstG})+(${cdsU}-${firstU}))." $cpd_map{ppi} ";

	    foreach my $factor (sort keys %{$factors{TranscriptionTerminationRNA}})
	    {					
		$rxn4 .= "+ 1 $factors{TranscriptionTerminationRNA}{$factor} ";
	    }
	    foreach my $factor (sort keys %{$factors{RNAP}})
	    {
		$rxn4 .= "+ 1 $factors{RNAP}{$factor} ";
	    }

	    $rxn4 .= "+ 3 $cpd_map{adp} + 3 $cpd_map{pi} + 3 $cpd_map{h}\tirreversible\tTranscription\t".((${cdsA}+${cdsC}+${cdsG}+${cdsU})/90)." (90nt/s) to ".((${cdsA}+${cdsC}+${cdsG}+${cdsU})/80)." (80nt/s) 1/s\n";
	    push @{$reactions{$gene}}, $rxn4;
	}
	else {
	    # print STDERR "skipping $gene of type $type\n";
	}

	# only translate CDS
	next unless ($type eq 'CDS');

	my ${cds_rna_C}=${cdsG}*10 + ${cdsC}*9 + ${cdsA}*10 + ${cdsU}*9;	
	my ${cds_rna_H}=${cdsG}*11 + ${cdsC}*11 + ${cdsA}*11 + ${cdsU}*10+1;
	my ${cds_rna_N}=${cdsG}*5 + ${cdsC}*3 + ${cdsA}*5 + ${cdsU}*2;
	my ${cds_rna_O}=${cdsG}*7 + ${cdsC}*7 + ${cdsA}*6 + ${cdsU}*8+7;
	my ${cds_rna_P}=(${cdsG} + ${cdsC}+${cdsU}+${cdsA})+2;
	my ${cds_rna_charge}=(${cdsG} + ${cdsC}+${cdsU}+${cdsA})*(-1)-3;

	push @{$reactions{$gene}}, "$gene\_mRNA_CONV\tconvsersion of mRNA to mRNA_1 (synthetic rxn)\t1 $gene\_mRNA --> 1 $gene\_mRNA_1\tirreversible\tTranslation\n";	
	
	my (%trnas_for_gene, $second_last_codon, $last_codon, $last_codon_);

	for (my $i=3; $i < (length($cds)-2); $i +=3) #code does not account for first triplet but set first aa always to be fmet;
	{
	    my $codon =substr($cds,$i,3);

	    if ($i < (length($cds)-3)) # MDJ: changed from -2 to -3 so we don't try to code an AA for the stop codon
	    {				
		if($fr =~ /programmed frameshift\-containing/ && $codon eq 'TGA')
		{
		    $i=$i+1;
		    $codon =substr($cds,$i,3);
		}

		if (exists $code_tRNA{$codon}{trna})
		{	
		    ++$trnas_for_gene{$code_tRNA{$codon}{trna}}; # counts the number of different type of tRNAs needed
		}
		else
		{
		    print STDERR "$gene\t$i\t$codon does not exist\n";
		}		
	    }

	    if ($i == (length($cds)-6))
	    {
		$second_last_codon=$code_tRNA{$codon}{trna};
	    }
	    elsif ($i == (length($cds)-3))
	    {
		$last_codon=$code_tRNA{$codon}{trna};
		$last_codon_=$codon; # keeps the last codon for each CDS - important for Translation (Release factor)
	    }
	}	

	if (! defined $last_codon && ! defined $second_last_codon) {
	    print STDERR "Skipping $gene because unable to identify last two codons\n";
	    next;
	}

	my $sum_aa = (length $aa_seq{$gene})-1; # subtracting one because ignoring first codon (fmet_tRNA or something like that)
	my ${fmet_tRNA}=1;
	my ${Mg2}=5*$sum_aa;
	my ${h2o}=2*$sum_aa-1;
	my ${gdp}=2*$sum_aa;
	my ${pi}=2*$sum_aa;
	my ${h}=2*$sum_aa+1;
	my $number_rib=int(($sum_aa+1)/17);
	
	my %aa_count = &aa_count($aa_seq{$gene}, \%infos);

	my $cpd_aa = "$gene\_aa\tpolypeptide $gene\tC$aa_count{C}H$aa_count{H}N$aa_count{N}O$aa_count{O}S$aa_count{S}";
	if ($aa_count{Se} !=0)
	{
	    $cpd_aa .= "Se$aa_count{Se}";
	}
	$cpd_aa .= "\t$aa_count{charge}\tTranslation\t\n";
	push @{$compounds{$gene}}, $cpd_aa;
	push @{$compounds{$gene}}, "$gene\_mRNA\tmRNA $gene\tC${cds_rna_C}H${cds_rna_H}N${cds_rna_N}O${cds_rna_O}P${cds_rna_P}\t${cds_rna_charge}\tTranslation\t".(${cdsG}+${cdsC}+${cdsU}+${cdsA})."(), , ${cdsG}G ${cdsC}C ${cdsU}U ${cdsA}A\n";
	push @{$compounds{$gene}}, "$gene\_mRNA_1\tmRNA $gene\tC${cds_rna_C}H${cds_rna_H}N${cds_rna_N}O${cds_rna_O}P${cds_rna_P}\t${cds_rna_charge}\tTranslation\t".(${cdsG}+${cdsC}+${cdsU}+${cdsA})."(), , ${cdsG}G ${cdsC}C ${cdsU}U ${cdsA}A\n";
	push @{$compounds{$gene}}, "$gene\_mRNA_2\tmRNA $gene\tC${cds_rna_C}H${cds_rna_H}N${cds_rna_N}O${cds_rna_O}P${cds_rna_P}\t${cds_rna_charge}\tTranslation\t".(${cdsG}+${cdsC}+${cdsU}+${cdsA})."(), , ${cdsG}G ${cdsC}C ${cdsU}U ${cdsA}A\n";
	
	#$gene\_mRNA_degr accounts for mRNA and degradosome (1 RNase E, 4 RhlB, 4 Eno, 4 PNPase) and Orn;
	my ${comp_degr_C}=${cds_rna_C};
	my ${comp_degr_H}=${cds_rna_H};
	my ${comp_degr_N}=${cds_rna_N};
	my ${comp_degr_O}=${cds_rna_O};
	my ${comp_degr_S}=0;
	my ${comp_degr_P}=${cds_rna_P};
	my ${comp_degr_Mg}=0;
	my ${comp_degr_Zn}=0;
	my ${comp_degr_Fe}=0;
	my ${comp_degr_charge}=${cds_rna_charge};
	my $NameFactors = '';			

	foreach my $factor (keys %{$factors{mRNAdegradation}})
	{				
	    ${comp_degr_C}+=$formulae{$factor}{C};
	    ${comp_degr_H}+=$formulae{$factor}{H};
	    ${comp_degr_N}+=$formulae{$factor}{N};
	    ${comp_degr_O}+=$formulae{$factor}{O};
	    ${comp_degr_S}+=$formulae{$factor}{S};
	    ${comp_degr_P}+=$formulae{$factor}{P};
	    ${comp_degr_Mg}+=$formulae{$factor}{Mg};
	    ${comp_degr_Zn}+=$formulae{$factor}{Zn};
	    ${comp_degr_Fe}+=$formulae{$factor}{Fe};
	    ${comp_degr_charge}+=$formulae{$factor}{charge};
	    $NameFactors.= $factor.", ";
	}
	
	push @{$compounds{$gene}}, "$gene\_mRNA_2_degr\tmRNA $gene degradation complex ($NameFactors)\tC${comp_degr_C}H${comp_degr_H}N${comp_degr_N}O${comp_degr_O}S${comp_degr_S}P${comp_degr_P}Mg${comp_degr_Mg}Zn${comp_degr_Zn}Fe${comp_degr_Fe}\t${comp_degr_charge}\tmRNA degradation\t, ${cdsG}G ${cdsC}C ${cdsU}U ${cdsA}A\n";

	my %as = map { $_ => 1 } (1, int($number_rib/2), $number_rib);

	foreach my $a (sort keys %as)
	{
	    ####################
	    # TRANSLATION INITIATION
	    ####################
	    my $rxn5 = "tl_ini_$gene\_$a\_rib\tTranslation initiation $gene $a ribosome(s) bound\t1 $gene\_mRNA_1 + $a fmet_tRNA_met ";
	    foreach my $factor (sort keys %{$factors{TranslationIni}}) {
		$rxn5 .= "+ $a $factors{TranslationIni}{$factor} ";
	    }				
	    $rxn5 .= "+ $a $cpd_map{h2o} --> 1 rib_ini_$gene\_$a ";
	    foreach my $factor (sort keys %{$factors{TranslationIniOut}})
	    {
		$rxn5 .= "+ $a $factors{TranslationIniOut}{$factor} ";
	    }		
	    $rxn5 .= "+ $a $cpd_map{pi} + $a $cpd_map{h}\tirreversible\tTranslation\n";
	    push @{$reactions{$gene}}, $rxn5;

	    #####################
	    #addition of formyl-methionine tRNA, rib_30 subunit
	    # COMPOSITION OF rib_ini_$gene\_$a
	    my $C_ini=${cds_rna_C}+$a*$formulae{fmet_tRNA_met}{C};
	    my $H_ini=${cds_rna_H}+$a*$formulae{fmet_tRNA_met}{H};
	    my $N_ini=${cds_rna_N}+$a*$formulae{fmet_tRNA_met}{N};
	    my $O_ini=${cds_rna_O}+$a*$formulae{fmet_tRNA_met}{O};
	    my $S_ini=$a*$formulae{fmet_tRNA_met}{S};
	    my $Se_ini=0;
	    my $Mg_ini=$a*$formulae{fmet_tRNA_met}{Mg};
	    my $Zn_ini=$a*$formulae{fmet_tRNA_met}{Zn};
	    my $Fe_ini=$a*$formulae{fmet_tRNA_met}{Fe};
	    my $P_ini=${cds_rna_P}+$a*$formulae{fmet_tRNA_met}{P};
	    my $Ch_ini=${cds_rna_charge}+$a*$formulae{fmet_tRNA_met}{charge};

	    my $NameFactors = '';			
	    foreach my $factor (sort keys %{$factors{TranslationIni}})
	    {
		$C_ini+=$a*$formulae{$factor}{C};
		$H_ini+=$a*$formulae{$factor}{H};
		$N_ini+=$a*$formulae{$factor}{N};
		$O_ini+=$a*$formulae{$factor}{O};
		$S_ini+=$a*$formulae{$factor}{S};
		$P_ini+=$a*$formulae{$factor}{P};
		$Mg_ini+=$a*$formulae{$factor}{Mg};
		$Zn_ini+=$a*$formulae{$factor}{Zn};
		$Fe_ini+=$a*$formulae{$factor}{Fe};
		$Se_ini+=$a*$formulae{$factor}{Se} if defined $formulae{$factor}{Se};
		$Ch_ini+=$a*$formulae{$factor}{charge};
		$NameFactors.= $factor.", ";
	    }
	    # substracts Initiation factors from initiation complex
	    foreach my $factor (sort keys %{$factors{TranslationIniFactors}})
	    {				
		$C_ini-=$a*$formulae{$factor}{C};
		$H_ini-=$a*$formulae{$factor}{H};
		$N_ini-=$a*$formulae{$factor}{N};
		$O_ini-=$a*$formulae{$factor}{O};
		$S_ini-=$a*$formulae{$factor}{S};
		$P_ini-=$a*$formulae{$factor}{P};
		$Mg_ini-=$a*$formulae{$factor}{Mg};
		$Zn_ini-=$a*$formulae{$factor}{Zn};
		$Fe_ini-=$a*$formulae{$factor}{Fe};
		$Se_ini-=$a*$formulae{$factor}{Se} if defined $formulae{$factor}{Se};
		$Ch_ini-=$a*$formulae{$factor}{charge};
	    }	
	    ####################
	    push @{$compounds{$gene}}, "rib_ini_$gene\_$a\t$gene complexed with $a* fmet-tRNA/Mg, $NameFactors\tC$C_ini"."H$H_ini"."N$N_ini"."O$O_ini"."S$S_ini"."Se$Se_ini"."Mg$Mg_ini"."P$P_ini"."Fe$Fe_ini"."Zn$Zn_ini\t$Ch_ini\tTranslation\n";
	    ####################

	    ####################
	    # TRANSLATION ELONGATION - PART Ia
	    ###################

	    my $rxn6 = "tl_elo_$gene\_$a\_rib1\tTranslation elongation 1 $gene $a ribosome(s)\t1 rib_ini_$gene\_$a + ".(${Mg2}*$a)." $cpd_map{mg2} ";

	    # preparation for elongation complex: rib_70_$gene\_$a\_cplx
	    # different elements come to the rib_ini complex : EF-TU.GTP.tRNA (charged with aa), EF-G(GTP)

	    my $C_elo=$C_ini;
	    my $H_elo=$H_ini;
	    my $N_elo=$N_ini;
	    my $O_elo=$O_ini;
	    my $S_elo=$S_ini;
	    my $Se_elo=$Se_ini;
	    my $Ch_elo=$Ch_ini+(${Mg2}*$a)*2;
	    my $P_elo=$P_ini;
	    my $Mg_elo=$Mg_ini+(${Mg2}*$a);
	    my $Zn_elo=$Zn_ini;
	    my $Fe_elo=$Fe_ini;

	    ###################
	    # PREPARATION FOR COMPLEX OF TRANSLATION ELONGATION PART II
	    ##################
	    # Accounts for mRNA + formylmethionine, ribosome 70S
	    my $C_elo2=${cds_rna_C}+$a*6;
	    my $H_elo2=${cds_rna_H}+$a*10;
	    my $N_elo2=${cds_rna_N}+$a*1;
	    my $O_elo2=${cds_rna_O}+$a*3;
	    my $P_elo2=${cds_rna_P}+$a*0;
	    my $S_elo2=0;
	    my $Se_elo2=0;
	    my $Mg_elo2=0;
	    my $Zn_elo2=0;
	    my $Fe_elo2=0;
	    my $Ch_elo2=${cds_rna_charge}+$a*(-1);
	    $NameFactors = '';			
	    foreach my $factor (sort keys %{$factors{TranslationElo2}})
	    {
		$C_elo2+=$a*$formulae{$factor}{C};
		$H_elo2+=$a*$formulae{$factor}{H};
		$N_elo2+=$a*$formulae{$factor}{N};
		$O_elo2+=$a*$formulae{$factor}{O};
		$S_elo2+=$a*$formulae{$factor}{S};
		$P_elo2+=$a*$formulae{$factor}{P};
		$Mg_elo2+=$a*$formulae{$factor}{Mg};
		$Zn_elo2+=$a*$formulae{$factor}{Zn};
		$Fe_elo2+=$a*$formulae{$factor}{Fe};
		$Se_elo2+=$a*$formulae{$factor}{Se} if defined $formulae{$factor}{Se};
		$Ch_elo2+=$a*$formulae{$factor}{charge};
		$NameFactors.= $factor.", ";
	    }				
	    
	    foreach my $trna (sort keys %trnas_for_gene)
	    {
		my $trna2;
		my $found = 0;
		# Can be only 1 factor
		foreach my $factor (sort keys %{$factors{TranslationEF_TU_GTP}})
		{
		    if ($found == 0)
		    {
			if ($trna =~ /^tRNA-(\w{3})-/) {
			    my @codons = keys %{$genetic_code{$assigned_code}{$1}};
			    $trna2 = $factor.'.'.$trna;
			    $found = 1;
			}
			else {
			    print STDERR "no genetic code match for $trna\n";
			}
		    }
		}

		################
		$rxn6 .= "+ ".($trnas_for_gene{$trna} * $a)." $trna2 ";
		################
		
		# addition of EF-Tu.GTP.aa-tRNA to initiation complex
		### DEFINE $trna2 formulae have to be predofined!
		$C_elo=$C_elo+(($trnas_for_gene{$trna})*$a*$formulae{$trna2}{C});
		$H_elo=$H_elo+(($trnas_for_gene{$trna})*$a*$formulae{$trna2}{H});
		$N_elo=$N_elo+(($trnas_for_gene{$trna})*$a*$formulae{$trna2}{N});
		$O_elo=$O_elo+(($trnas_for_gene{$trna})*$a*$formulae{$trna2}{O});
		$S_elo=$S_elo+(($trnas_for_gene{$trna})*$a*$formulae{$trna2}{S});
		$Se_elo=$Se_elo+(($trnas_for_gene{$trna})*$a*$formulae{$trna2}{Se}) if defined $formulae{$trna2}{Se};
		$Ch_elo=$Ch_elo+(($trnas_for_gene{$trna})*$a*$formulae{$trna2}{charge});
		$P_elo=$P_elo+($trnas_for_gene{$trna})*$a*$formulae{$trna2}{P};
		$Mg_elo=$Mg_elo+($trnas_for_gene{$trna})*$a*$formulae{$trna2}{Mg};
		$Zn_elo=$Zn_elo+($trnas_for_gene{$trna})*$a*$formulae{$trna2}{Zn};
		$Fe_elo=$Fe_elo+($trnas_for_gene{$trna})*$a*$formulae{$trna2}{Fe};
		
		##################
		# accounts for fact that $a*last tRNA and $a*EF-TU is still attached to complex and $a*EF-G
		##################
		if ($trna eq ${second_last_codon})
		{
		    $C_elo2=$C_elo2+$infos{$list_tRNA{$trna}}{C}*(($trnas_for_gene{$trna}-1) * $a)+($a*$formulae{$trna2}{C});
		    $H_elo2=$H_elo2+$infos{$list_tRNA{$trna}}{H}*(($trnas_for_gene{$trna}-1) * $a)+($a*$formulae{$trna2}{H});
		    $N_elo2=$N_elo2+$infos{$list_tRNA{$trna}}{N}*(($trnas_for_gene{$trna}-1) * $a)+($a*$formulae{$trna2}{N});
		    $O_elo2=$O_elo2+$infos{$list_tRNA{$trna}}{O}*(($trnas_for_gene{$trna}-1) * $a)+($a*$formulae{$trna2}{O});
		    $P_elo2=$P_elo2+$infos{$list_tRNA{$trna}}{P}*(($trnas_for_gene{$trna}-1) * $a)+($a*$formulae{$trna2}{P}) if defined $infos{$list_tRNA{$trna}}{P};
		    $S_elo2=$S_elo2+$infos{$list_tRNA{$trna}}{S}*(($trnas_for_gene{$trna}-1) * $a)+($a*$formulae{$trna2}{S}) if defined $infos{$list_tRNA{$trna}}{S};
		    $Se_elo2=$Se_elo2+$infos{$list_tRNA{$trna}}{Se}*(($trnas_for_gene{$trna}-1) * $a)+($a*$formulae{$trna2}{Se}) if defined $infos{$list_tRNA{$trna}}{Se};
		    $Ch_elo2=$Ch_elo2+$infos{$list_tRNA{$trna}}{charge}*(($trnas_for_gene{$trna}-1) * $a)+($a*$formulae{$trna2}{charge});
		    $Mg_elo2=$Mg_elo2+$infos{$list_tRNA{$trna}}{Mg}*(($trnas_for_gene{$trna}-1) * $a)+($a*$formulae{$trna2}{Mg}) if defined $infos{$list_tRNA{$trna}}{Mg};
		    $Zn_elo2=$Zn_elo2+$infos{$list_tRNA{$trna}}{Zn}*(($trnas_for_gene{$trna}-1) * $a)+($a*$formulae{$trna2}{Zn}) if defined $infos{$list_tRNA{$trna}}{Zn};
		    $Fe_elo2=$Fe_elo2+$infos{$list_tRNA{$trna}}{Fe}*(($trnas_for_gene{$trna}-1) * $a)+($a*$formulae{$trna2}{Fe}) if defined $infos{$list_tRNA{$trna}}{Fe};
		    $NameFactors = '';			
		    foreach my $factor (sort keys %{$factors{TranslationEF_G_GTP}})
		    {
			$C_elo2+=$a*$formulae{$factor}{C};
			$H_elo2+=$a*$formulae{$factor}{H};
			$N_elo2+=$a*$formulae{$factor}{N};
			$O_elo2+=$a*$formulae{$factor}{O};
			$S_elo2+=$a*$formulae{$factor}{S};
			$P_elo2+=$a*$formulae{$factor}{P};
			$Mg_elo2+=$a*$formulae{$factor}{Mg};
			$Zn_elo2+=$a*$formulae{$factor}{Zn};
			$Fe_elo2+=$a*$formulae{$factor}{Fe};
			$Se_elo2+=$a*$formulae{$factor}{Se} if defined $formulae{$factor}{Se};
			$Ch_elo2+=$a*$formulae{$factor}{charge};
			$NameFactors.= $factor.", ";
		    }	
		}
		else
		{				
		    # PREPARATION FOR COMPLEX OF TRANSLATION ELONGATION PART II - adds aa to formulae
		    $C_elo2=$C_elo2+$infos{$list_tRNA{$trna}}{C}*($trnas_for_gene{$trna} * $a);
		    $H_elo2=$H_elo2+$infos{$list_tRNA{$trna}}{H}*($trnas_for_gene{$trna} * $a);
		    $N_elo2=$N_elo2+$infos{$list_tRNA{$trna}}{N}*($trnas_for_gene{$trna} * $a);
		    $O_elo2=$O_elo2+$infos{$list_tRNA{$trna}}{O}*($trnas_for_gene{$trna} * $a);
		    $S_elo2=$S_elo2+$infos{$list_tRNA{$trna}}{S}*($trnas_for_gene{$trna} * $a) if defined $infos{$list_tRNA{$trna}}{S};
		    $Se_elo2=$Se_elo2+$infos{$list_tRNA{$trna}}{Se}*($trnas_for_gene{$trna} * $a) if defined $infos{$list_tRNA{$trna}}{Se};
		    $Mg_elo2=$Mg_elo2+$infos{$list_tRNA{$trna}}{Mg}*($trnas_for_gene{$trna} * $a) if defined $infos{$list_tRNA{$trna}}{Mg};
		    $Ch_elo2=$Ch_elo2+$infos{$list_tRNA{$trna}}{charge}*($trnas_for_gene{$trna} * $a);
		}
	    }
	    
	    $NameFactors = '';			
	    foreach my $factor (sort keys %{$factors{TranslationEF_G_GTP}})
	    {
		$C_elo+=(${sum_aa}*$a)*$formulae{$factor}{C};
		$H_elo+=(${sum_aa}*$a)*$formulae{$factor}{H};
		$N_elo+=(${sum_aa}*$a)*$formulae{$factor}{N};
		$O_elo+=(${sum_aa}*$a)*$formulae{$factor}{O};
		$S_elo+=(${sum_aa}*$a)*$formulae{$factor}{S};
		$P_elo+=(${sum_aa}*$a)*$formulae{$factor}{P};
		$Mg_elo+=(${sum_aa}*$a)*$formulae{$factor}{Mg};
		$Zn_elo+=(${sum_aa}*$a)*$formulae{$factor}{Zn};
		$Fe_elo+=(${sum_aa}*$a)*$formulae{$factor}{Fe};
		$Se_elo+=(${sum_aa}*$a)*$formulae{$factor}{Se} if defined $formulae{$factor}{Se};
		$Ch_elo+=(${sum_aa}*$a)*$formulae{$factor}{charge};
		$NameFactors.= $factor.", ";
	    }		
	    
	    #################
	    # TRANSLATION ELONGATION RXN - PART Ib
	    #################
	    foreach my $factor (sort keys %{$factors{TranslationEF_G_GTP}})
	    {
		$rxn6 .= "+ ".(${sum_aa}*$a)." $factors{TranslationEF_G_GTP}{$factor} ";
	    }
	    $rxn6 .= "--> 1 rib_70_elo1_$gene\_$a\_cplx\tirreversible\tTranslation\n";
	    push @{$reactions{$gene}}, $rxn6;

	    #################
	    # WRITING OF rib_70_$gene\_$a\_cplx COMPOUND
	    #################

	    my $cpd_rib_70 = "rib_70_elo1_$gene\_$a\_cplx\tTranslation elongation complex: $a * ribosome 70S/$gene/EF-TU-tRNA's/EF-G\tC".$C_elo."H".$H_elo."N".$N_elo."O".$O_elo."S".$S_elo."";
	    if ($Se_elo!=0)
	    {
		$cpd_rib_70 .= "Se$Se_elo";
	    }
	    $cpd_rib_70 .= "P".$P_elo."Mg".$Mg_elo."Zn".$Zn_elo."Fe".$Fe_elo."\t$Ch_elo\tTranslation\n";
	    push @{$compounds{$gene}}, $cpd_rib_70;
	    
	    #################
	    # WRITING OF rib_70_elo2_$gene\_$a\_cplx COMPOUND
	    #################

	    # substracts H2O from resulting complex since H2O is released when peptide bond is formed
	    $H_elo2=$H_elo2 - (${sum_aa}-1) * $a * 2;
	    $O_elo2=$O_elo2 - (${sum_aa}-1) * $a *1;

	    my $cpd_rib_70_2 = "rib_70_elo2_$gene\_$a\_cplx\tTranslation elongation complex: $a * ribosome 70S/$gene/second last codon EF-TU-tRNA/1 EF-G\tC".$C_elo2."H".$H_elo2."N".$N_elo2."O".$O_elo2."S".$S_elo2."";
	    if ($Se_elo2!=0)
	    {
		$cpd_rib_70_2 .= "Se$Se_elo2";
	    }
	    $cpd_rib_70_2 .= "P".$P_elo2."Mg".$Mg_elo2."Zn".$Zn_elo2."Fe".$Fe_elo2."\t$Ch_elo2\tTranslation\n";
	    push @{$compounds{$gene}}, $cpd_rib_70_2;

	    ##################
	    # TRANSLATION ELONGATION - PART II
	    ##################
	    
	    ####
	    # WRITE reaciton tl_elo_$gene\_$a\_rib2
	    # needs $factors{$factor}{TranslationEF-TU.GDP}
	    # needs $factors{$factor}{TranslationEF-G.GDP}
	    ####
	    ### THINK 
	    my $rxn7 = "tl_elo_$gene\_$a\_rib2\tTranslation elongation 2 $gene $a ribosome(s)\t1 rib_70_elo1_$gene\_$a\_cplx + ".(((${sum_aa}-1)*$a))." $cpd_map{h2o} --> 1 rib_70_elo2_$gene\_$a\_cplx + $a fmet_tRNA + ".(${Mg2}*$a)." $cpd_map{mg2} + ".(((2*${sum_aa})-2)*$a)." $cpd_map{pi} + ".(((2*${sum_aa})-2)*$a)." $cpd_map{h} ";
	    foreach my $factor (sort keys %{$factors{TranslationEF_TU_GDP}})
	    {
		$rxn7 .= "+ ".((${sum_aa}*$a)-$a)." $factors{TranslationEF_TU_GDP}{$factor} ";
	    }
	    foreach my $factor (sort keys %{$factors{TranslationEF_G_GDP}})
	    {
		$rxn7 .= "+ ".((${sum_aa}*$a)-$a)." $factors{TranslationEF_G_GDP}{$factor} ";
	    }
	    
	    foreach my $trna (keys %list_tRNA)
	    {			
		if (defined $trnas_for_gene{$trna} && $trnas_for_gene{$trna} != 0)
		{
		    if ($trna eq ${second_last_codon})
		    {		
			$rxn7 .= "+ ".($a*$trnas_for_gene{$trna}-$a); # does not remove last aa-trna from complex
			$rxn7 .= " $trna ";
		    }
		    else
		    {
			$rxn7 .= "+ ".($a*$trnas_for_gene{$trna});
			$rxn7 .= " $trna ";
		    }
		}
	    }
	    $rxn7 .= "\tirreversible\tTranslation\n";
	    push @{$reactions{$gene}}, $rxn7;
	    
	    #####################
	    # Translation TERMINATION
	    #####################
	    
	    # accounting for different release factors depending on last codon (based on Putzer abd Laalami, 2003, book, Translation Mechanisms,chapter 24)
	    my $rfx ='';
	    my $rf = '';
	    my $CountRF = 0;

	    foreach my $factor (sort keys %{$factors{LastCodonsFactors}})
	    {
		if (exists $factors{LastCodonsFactors}{$factor}{$last_codon_})
		{
		    if ($CountRF == 0)
		    {
			$rf=$factors{LastCodonsFactors}{$factor}{$last_codon_};
			++$CountRF;
		    }
		    else
		    {
			$rfx=$factors{LastCodonsFactors}{$factor}{$last_codon_};
		    }
		}
	    }
	    
	    ####
	    # Write tl_term_$gene\_$a\_rib1
	    # needs  $factors{RF3_mono.GDP}{TranslationTerm} = 1;
	    # $a RF3_mono.GDP + $a Rrf_mono ";
	    ####
	    
	    my $rxn8 = "tl_term_$gene\_$a\_rib1\tTranslation termination 1 $gene $a ribosome(s), ($rf)\t1 rib_70_elo2_$gene\_$a\_cplx ";
	    
	    foreach my $factor (sort keys %{$factors{TranslationTerm}})
	    {
		$rxn8 .= "+ $a $factors{TranslationTerm}{$factor} ";
	    }				
	    #THINK 
	    $rxn8 .= "+ $a $rf --> 1 tl_term_$gene\_$a\_rib_cplx\treversible\tTranslation\n";
	    push @{$reactions{$gene}}, $rxn8;

	    ####
	    # Write tl_term_$gene\_$a\_rib2
	    # needs $factors{rib_70}factors{TranslationElo2} = 1;
	    #+ $a rib_70 
	    # needs $factors{EF-Tu.GDP}{TranslationEF-TU.GDP} = 1;
	    #+ $a EF-Tu.GDP
	    # needs $factors{EF-G.GDP}{TranslationEF-G.GDP} = 1;
	    #+ $a EF-G.GDP
	    # needs  $factors{RF3_mono.GDP}{TranslationTerm} = 1;
	    # $a RF3_mono.GDP + $a Rrf_mono ";
	    ####				
	    
	    my $rxn9 = "tl_term_$gene\_$a\_rib2\tTranslation termination 2 $gene $a ribosome(s), ($rf)\t1 tl_term_$gene\_$a\_rib_cplx + $a $cpd_map{gtp} + ".(2*$a)." $cpd_map{h2o} --> 1 $gene\_mRNA_2 + $a $gene\_aa + $a $cpd_map{gdp} + ".(3*$a)." $cpd_map{h} + ".(3*$a)." $cpd_map{pi} + $a $rf ";
	    #+ $a rib_70 
	    foreach my $factor (sort keys %{$factors{TranslationElo2}})
	    {
		$rxn9 .= "+ $a $factors{TranslationElo2}{$factor} ";
	    }				
	    #+ $a EF-Tu.GDP
	    foreach my $factor (sort keys %{$factors{TranslationEF_TU_GDP}})
	    {
		$rxn9 .= "+ $a $factors{TranslationEF_TU_GDP}{$factor} ";
	    }				
	    #+ $a EF-G.GDP
	    foreach my $factor (sort keys %{$factors{TranslationEF_G_GDP}})
	    {
		$rxn9 .= "+ $a $factors{TranslationEF_G_GDP}{$factor} ";
	    }				
	    #+ $a RF3_mono.GDP + $a Rrf_mono  
	    foreach my $factor (sort keys %{$factors{TranslationTerm}})
	    {
		$rxn9 .= "+ $a $factors{TranslationTerm}{$factor} ";
	    }
	    $rxn9 .= "+ $a ";
	    if (defined ${second_last_codon})	
	    {
		my $trna=${second_last_codon};
		$rxn9 .= "$trna\tirreversible\tTranslation\n";
		push @{$reactions{$gene}}, $rxn9;
	    }

	    if ($rfx ne '')
	    {
		my $rxn11 = "tl_term_$gene\_$a\_rib_RF2_1\tTranslation termination 1 $gene $a ribosome(s), RF2_mono\t1 rib_70_elo2_$gene\_$a\_cplx + $a $rfx ";
		#+ $a RF3_mono.GDP + $a Rrf_mono 
		foreach my $factor (sort keys %{$factors{TranslationTerm}})
		{
		    $rxn11 .= "+ $a $factors{TranslationTerm}{$factor} ";
		}					
		$rxn11 .= "--> 1 tl_term_$gene\_$a\_rib_cplx2\treversible\tTranslation\n";
		push @{$reactions{$gene}}, $rxn11;
		
		my $rxn22 = "tl_term_$gene\_$a\_rib_RF2_2\tTranslation termination 2 $gene $a ribosome(s), RF2\t1 tl_term_$gene\_$a\_rib_cplx2 + $a $cpd_map{gtp} + ".(2*$a)." $cpd_map{h2o} --> 1 $gene\_mRNA_2 + $a $gene\_aa + $a $cpd_map{gdp} + ".(3*$a)." $cpd_map{h} + ".(3*$a)." $cpd_map{pi} + $a $rfx ";

		#+ $a rib_70 
		foreach my $factor (sort keys %{$factors{TranslationElo2}})
		{
		    $rxn22 .= "+ $a $factors{TranslationElo2}{$factor} ";
		}				
		#+ $a EF-Tu.GDP
		foreach my $factor (sort keys %{$factors{TranslationEF_TU_GDP}})
		{
		    $rxn22 .= "+ $a $factors{TranslationEF_TU_GDP}{$factor} ";
		}				
		#+ $a EF-G.GDP
		foreach my $factor (sort keys %{$factors{TranslationEF_G_GDP}})
		{
		    $rxn22 .= "+ $a $factors{TranslationEF_G_GDP}{$factor} ";
		}				
		#+ $a RF3_mono.GDP + $a Rrf_mono  
		foreach my $factor (sort keys %{$factors{TranslationTerm}})
		{
		    $rxn22 .= "+ $a $factors{TranslationTerm}{$factor} ";
		}					
		$rxn22 .= "+ $a ";
		my $trna=${second_last_codon};
		$rxn22 .= "$trna\tirreversible\tTranslation\n";
		push @{$reactions{$gene}}, $rxn22;
		
		###################
		# CREATES tl_term_$gene\_$a\_rib_cplx complex, contains: rib_70_elo2_$gene\_$a, $a*RF3_mono.GDP, $a*Rrf_mono, $a*rfx (RF2_mono)
		###################

		my $C_term2=$C_elo2 + $a*$formulae{$rfx}{C};
		my $H_term2=$H_elo2 + $a*$formulae{$rfx}{H};
		my $N_term2=$N_elo2 + $a*$formulae{$rfx}{N};
		my $O_term2=$O_elo2 + $a*$formulae{$rfx}{O};
		my $P_term2=$P_elo2 + $a*$formulae{$rfx}{P};
		my $S_term2=$S_elo2 + $a*$formulae{$rfx}{S};
		my $Se_term2=$Se_elo2; $Se_term2=$Se_term2 + $a*$formulae{$rfx}{Se} if defined $formulae{$rfx}{Se};
		my $Ch_term2=$Ch_elo2 + $a*$formulae{$rfx}{charge};
		my $Mg_term2=$Mg_elo2 + $a*$formulae{$rfx}{Mg};
		my $Zn_term2=$Zn_elo2 + $a*$formulae{$rfx}{Zn};
		my $Fe_term2=$Fe_elo2 + $a*$formulae{$rfx}{Fe};
		
		$NameFactors = '';			
		foreach my $factor (sort keys %{$factors{TranslationTerm}})
		{
		    $C_term2+=$a*$formulae{$factor}{C};
		    $H_term2+=$a*$formulae{$factor}{H};
		    $N_term2+=$a*$formulae{$factor}{N};
		    $O_term2+=$a*$formulae{$factor}{O};
		    $S_term2+=$a*$formulae{$factor}{S};
		    $P_term2+=$a*$formulae{$factor}{P};
		    $Mg_term2+=$a*$formulae{$factor}{Mg};
		    $Zn_term2+=$a*$formulae{$factor}{Zn};
		    $Fe_term2+=$a*$formulae{$factor}{Fe};
		    $Se_term2+=$a*$formulae{$factor}{Se} if defined $formulae{$factor}{Se};
		    $Ch_term2+=$a*$formulae{$factor}{charge};
		    $NameFactors.= $factor.", ";
		}

		my $cpd_cplx2 = "tl_term_$gene\_$a\_rib_cplx2\tTranslation termination complex $gene $a ribosomes 2 ($rfx)\tC".$C_term2."H".$H_term2."N".$N_term2."O".$O_term2."S".$S_term2."";
		if($Se_term2!=0)
		{
		    $cpd_cplx2 .= "Se$Se_term2";
		}
		$cpd_cplx2 .= "P".$P_term2."Mg".$Mg_term2."Zn".$Zn_term2."Fe".$Fe_term2."\t$Ch_term2\tTranslation\n";
		push @{$compounds{$gene}}, $cpd_cplx2;
	    }

	    ###################
	    # CREATES tl_term_$gene\_$a\_rib_cplx complex, contains: rib_70_elo2_$gene\_$a, $a*RF3_mono.GDP, $a*Rrf_mono, $a*rf (RF1 or RF2)
	    ###################
	    my $C_term1=$C_elo2 + $a*$formulae{$rf}{C};
	    my $H_term1=$H_elo2 + $a*$formulae{$rf}{H};
	    my $N_term1=$N_elo2 + $a*$formulae{$rf}{N};
	    my $O_term1=$O_elo2 + $a*$formulae{$rf}{O};
	    my $P_term1=$P_elo2 + $a*$formulae{$rf}{P};
	    my $S_term1=$S_elo2 + $a*$formulae{$rf}{S};
	    my $Se_term1=$Se_elo2; $Se_term1 += $a*$formulae{$rf}{Se} if defined $formulae{$rf}{Se};
	    my $Ch_term1=$Ch_elo2 + $a*$formulae{$rf}{charge};
	    my $Mg_term1=$Mg_elo2 + $a*$formulae{$rf}{Mg};
	    my $Zn_term1=$Zn_elo2 + $a*$formulae{$rf}{Zn};
	    my $Fe_term1=$Fe_elo2 + $a*$formulae{$rf}{Fe};
	    
	    $NameFactors = '';			
	    foreach my $factor (sort keys %{$factors{TranslationTerm}})
	    {
		$C_term1+=$a*$formulae{$factor}{C};
		$H_term1+=$a*$formulae{$factor}{H};
		$N_term1+=$a*$formulae{$factor}{N};
		$O_term1+=$a*$formulae{$factor}{O};
		$S_term1+=$a*$formulae{$factor}{S};
		$P_term1+=$a*$formulae{$factor}{P};
		$Mg_term1+=$a*$formulae{$factor}{Mg};
		$Zn_term1+=$a*$formulae{$factor}{Zn};
		$Fe_term1+=$a*$formulae{$factor}{Fe};
		$Se_term1+=$a*$formulae{$factor}{Se} if defined $formulae{$factor}{Se};
		$Ch_term1+=$a*$formulae{$factor}{charge};
		$NameFactors.= $factor.", ";
	    }

	    my $cpd_term = "tl_term_$gene\_$a\_rib_cplx\tTranslation termination complex $gene $a ribosomes ($NameFactors)\tC".$C_term1."H".$H_term1."N".$N_term1."O".$O_term1."S".$S_term1."";
	    if($Se_term1!=0)
	    {
		$cpd_term .= "Se$Se_term1";
	    }
	    $cpd_term .= "P".$P_term1."Mg".$Mg_term1."Zn".$Zn_term1."Fe".$Fe_term1."\t$Ch_term1\tTranslation\n";
	    push @{$compounds{$gene}}, $cpd_term;
	}

	# mRNA degradation
	my $cds_length=int(($cdsA+$cdsC+$cdsG+$cdsU)*0.25); 
	
	###
	# WRITE $cds\_mRNA_degr1
	# needs $factors{degradosome}{mRNAdegradation} = 0;
	#+ 1 degradosome + 1 Orn_dim 
	###
	my $rxn33 = "$gene\_mRNA_degr1\tDegradation of $gene\_mRNA I\t1 $gene\_mRNA_2 ";

	foreach my $factor (sort keys %{$factors{mRNAdegradation}})
	{				
	    $rxn33 .= "+ 1 ".($factors{mRNAdegradation}{$factor})." ";
	}		
	$rxn33 .= "+ $cds_length $cpd_map{atp} + $cds_length $cpd_map{h2o} --> 1 $gene\_mRNA_2_degr + $cds_length $cpd_map{adp} + $cds_length $cpd_map{pi} + $cds_length $cpd_map{h} \treversible\tmRNA degradation\n"; #accounts for unwinding of ds mRNA; based on Gralla and DiLisi, 1974, Nature, 248,330-332 --> found that 50% of randomly generated mRNA has to be ds mRNA
	push @{$reactions{$gene}}, $rxn33;
	
	###
	# WRITE $cds\_mRNA_degr2
	# needs $factors{degradosome}{mRNAdegradation} = 0;
	#+ 1 degradosome + 1 Orn_dim 
	###		
	
	my $rxn44 = "$gene\_mRNA_degr2\tDegradation of $gene\_mRNA II\t1 $gene\_mRNA_2_degr + ".(${cdsA}+${cdsC}+${cdsG}+${cdsU}-1)." $cpd_map{h2o} --> ";
	my $cnt = 0;
	foreach my $factor (sort keys %{$factors{mRNAdegradation}})
	{		
	    if ($cnt == 0)
	    {
		$rxn44 .= "1 ".($factors{mRNAdegradation}{$factor})." ";
		$cnt +=1;
	    }
	    else
	    {
		$rxn44 .= "+ 1 ".($factors{mRNAdegradation}{$factor})." ";
	    }
	}
	
	if (substr($cds,0,1) eq 'A')
	{
	    $rxn44 .= "+ 1 ".$cpd_map{atp}." + ".(${cdsA}-1)." ".$cpd_map{amp};
	    $rxn44 .= " + ${cdsC} ".$cpd_map{cmp}." + ${cdsG} ".$cpd_map{gmp}." + ${cdsU} ".$cpd_map{ump}." + ".(${cdsA}+${cdsC}+${cdsG}+${cdsU}-1)." ".$cpd_map{h}."\tirreversible\tmRNA degradation\n";
	}
	elsif (substr($cds,0,1) eq 'G')
	{
	    $rxn44 .= "+ 1 ".$cpd_map{gtp}." + ".(${cdsG}-1)." ".$cpd_map{gmp};
	    $rxn44 .= " + ${cdsC} ".$cpd_map{cmp}." + ${cdsA} ".$cpd_map{amp}." + ${cdsU} ".$cpd_map{ump}." + ".(${cdsA}+${cdsC}+${cdsG}+${cdsU}-1)." ".$cpd_map{h}."\tirreversible\tmRNA degradation\n";
	}
	elsif(substr($cds,0,1) eq 'T')
	{
	    $rxn44 .= "+ 1 ".$cpd_map{utp}." + ".(${cdsU}-1)." ".$cpd_map{ump};
	    $rxn44 .= " + ${cdsC} ".$cpd_map{cmp}." + ${cdsG} ".$cpd_map{gmp}." + ${cdsA} ".$cpd_map{amp}." + ".(${cdsA}+${cdsC}+${cdsG}+${cdsU}-1)." ".$cpd_map{h}."\tirreversible\tmRNA degradation\n";
	}
	elsif(substr($cds,0,1) eq 'C')
	{
	    $rxn44 .= "+ 1 ".$cpd_map{ctp}." + ".(${cdsC}-1)." ".$cpd_map{cmp};
	    $rxn44 .= " + ${cdsU} ".$cpd_map{ump}." + ${cdsG} ".$cpd_map{gmp}." + ${cdsA} ".$cpd_map{amp}." + ".(${cdsA}+${cdsC}+${cdsG}+${cdsU}-1)." ".$cpd_map{h}."\tirreversible\tmRNA degradation\n";
	}
	else
	{
	    $rxn44 .= "\tmRNA degradation\n";
	}
	push @{$reactions{$gene}}, $rxn44;

	###
	# WRITE $gen\_maturation1
	#  needs $factors{Def_mono}{ProtMaturationDef} = 1;
	#+ 1 Def_mono 
	###	
	my $rxn55 = "$gene\_maturation1\tPolypeptide $gene peptide deformylase complex\t1 $gene\_aa ";
	foreach my $factor (keys %{$factors{ProtMaturationDef}})
	{				
	    $rxn55 .= "+ 1 $factors{ProtMaturationDef}{$factor} ";
	}
	$rxn55 .= "--> 1 $gene\_def_cplx\treversible\tProtein Maturation\tno matured protein available-missing entry\n";
	push @{$reactions{$gene}}, $rxn55;

	###
	# WRITE $gene\_maturation2
	#  needs $factors{Def_mono}{ProtMaturationDef} = 1;
	#+ 1 Def_mono 
	###		
	my $rxn66 = "$gene\_maturation2\t$gene formation\t1 $gene\_def_cplx + 1 $cpd_map{h2o} --> 1 $gene\_m ";
	foreach my $factor (keys %{$factors{ProtMaturationDef}})
	{				
	    $rxn66 .= "+ 1 $factors{ProtMaturationDef}{$factor} ";
	}
	$rxn66 .= "+ 1 $cpd_map{for}\tirreversible\tProtein Maturation\tno matured protein available-missing entry\n";
	push @{$reactions{$gene}}, $rxn66;
	
	my ${mat_C} = $aa_count{C}; 			
	my ${mat_H} = $aa_count{H}; 
	my ${mat_O} = $aa_count{O}; 
	my ${mat_N} = $aa_count{N};
	my ${mat_S} = $aa_count{S}; 
	my ${mat_Mg} = 0;
	my ${mat_Zn} = 0;
	my ${mat_Fe} = 0;
	my ${mat_Se} = $aa_count{Se}; 
	my ${mat_charge} = $aa_count{charge}; 
	
	foreach my $factor (keys %{$factors{ProtMaturationDef}})
	{				
	    ${mat_C} += $formulae{$factor}{C}; 
	    ${mat_H} += $formulae{$factor}{H}; 
	    ${mat_N} += $formulae{$factor}{N}; 
	    ${mat_O} += $formulae{$factor}{O}; 
	    ${mat_S} += $formulae{$factor}{S}; 
	    ${mat_Mg} += $formulae{$factor}{Mg}; 
	    ${mat_Zn} += $formulae{$factor}{Zn}; 
	    ${mat_Fe} += $formulae{$factor}{Fe}; 
	    ${mat_Se} += $formulae{$factor}{Se} if exists $formulae{$factor}{Se};
	    ${mat_charge} += $formulae{$factor}{charge}; 
	}

	my $cpd_def = "$gene\_def_cplx\tPolypeptide $gene peptide deformylase complex\tC${mat_C}H${mat_H}N${mat_N}O${mat_O}S${mat_S}";
	if ($aa_count{Se}!=0)
	{
	    $cpd_def .= "Se${mat_Se}";
	}
	$cpd_def .= "Mg${mat_Mg}Zn${mat_Zn}Fe${mat_Fe}\t${mat_charge}\tMaturation\t\n";
	push @{$compounds{$gene}}, $cpd_def;

	#### sets bxxxx_m to bxxxx_aa
	my %m_count = &mature_count(%aa_count);
	my $mw_protein=$mw{C}*$m_count{C}+$mw{H}*$m_count{H}+$mw{N}*$m_count{N}+$mw{O}*$m_count{O}+$mw{S}*$m_count{S};
	push @{$compounds{$gene}}, "$gene\_m\tMatured polypeptide $gene\tC".($m_count{C})."H".($m_count{H})."N$m_count{N}O".($m_count{O})."S$m_count{S}\t".($m_count{charge})."\tMaturation\tMW: $mw_protein\n";

	push @{$reactions{$gene}}, "$gene\_fold_spon\t$gene\_m folding: spontanous\t1 $gene\_m --> 1 $gene\_mono\tirreversible\tProtein Folding\n";    

	$mw_protein=$mw{C}*$m_count{C}+$mw{H}*$m_count{H}+$mw{N}*$m_count{N}+$mw{O}*$m_count{O}+$mw{S}*$m_count{S};
	push @{$compounds{$gene}}, "$gene\_mono\tMonomer\tC".($m_count{C})."H".($m_count{H})."N$m_count{N}O".($m_count{O})."S$m_count{S}\t".($m_count{charge})."\tFolding\tMW: $mw_protein\n";
    }

# MDJ: Add recycling reactions
    foreach my $factor (keys %{$factors{TranslationEF_TU_GDP}}) {
	# leave EF-Tu-Ts because it is cycled intermediate
	push @{$reactions{"Recycling"}}, "EF_Tu_cycle_1\tEF-Tu.GDP dissociation with EF-Ts as intermediary\t1 $factors{TranslationEF_TU_GDP}{$factor} + 1 $factors{TranslationEF_TS}{'EF-Ts'} --> 1 EF-Tu-Ts + 1 $cpd_map{gdp}\tirreversible\tRecycling\n";
    }
    foreach my $factor (keys %{$factors{TranslationEF_TU_GTP}}) {
	push @{$reactions{"Recycling"}}, "EF_Tu_cycle_2\tEF-Tu-Ts dissociation with GTP charging\t1 EF-Tu-Ts + 1 $cpd_map{gtp} --> 1 $factors{TranslationEF_TU_GTP}{$factor} + 1 $factors{TranslationEF_TS}{'EF-Ts'}\tirreversible\tRecycling\n";
    }

    foreach my $trna (keys %list_tRNA) {
	if ($trna =~ /^tRNA\-(\w{3}(\(p\))?)\-\w{3}/) {
	    my $aa = $1;
	    my @codons = keys %{$genetic_code{$assigned_code}{$aa}};
	    next if @codons == 0;
	    my $trna2 = $genetic_code{$assigned_code}{$1}{$codons[0]};
	    foreach my $factor (sort keys %{$factors{TranslationEF_TU_GTP}}) {
		push @{$reactions{"Recycling"}}, "EF_Tu_cycle_3_$trna\tCharging EF-Tu with $trna and GTP and AA\t1 $factors{TranslationEF_TU_GTP}{$factor} + 1 $trna + 1 $cpd_map{$aa} --> 1 EF-Tu.GTP.$trna\tirreversible\tRecycling\n";
	    }
	}
    }

    foreach my $factor (sort keys %{$factors{TranslationEF_G_GDP}})
    {
	foreach my $factor2 (sort keys %{$factors{TranslationEF_G_GTP}})
	{
	    push @{$reactions{"Recycling"}}, "EF_G_cycle\tEF-G.GDP dissociation and reassociation with GTP\t1 $factors{TranslationEF_G_GDP}{$factor} + 1 $cpd_map{gtp} --> 1 $factors{TranslationEF_G_GTP}{$factor2} + 1 $cpd_map{gdp}\tirreversible\tRecycling\n";
	}				
    }				

# not sure how this happens exactly, but probably OK - NEED TO MASS BALANCE AND ADD COFACTORS (E.G., H2O?)
# leave IF2 all by itself because it only appears in recycling reactions
    push @{$reactions{"Recycling"}}, "IF2.GDP_cycle\tIF2.GDP dissociation\t1 $factors{TranslationIniOut}{'IF2.GDP'} --> 1 IF2 + 1 $cpd_map{gdp}\tirreversible\tRecycling\n";
    push @{$reactions{"Recycling"}}, "Rib_cycle\tInactive 70S ribosome dissociation plus initiation factor binding\t1 $factors{TranslationElo2}{rib_70} + 1 $factors{TranslationIniFactors}{IF1} + 1 IF2 + 1 $factors{TranslationIniFactors}{IF3} + $cpd_map{gtp} --> 1 $factors{TranslationIni}{rib_30_ini} + 1 $factors{TranslationIni}{rib_50}\tirreversible\tRecycling\n";
    foreach my $factor (keys %{$factors{RNAP}}) {
	push @{$reactions{"Recycling"}}, "RNAP_sigma_cycle_$factor\tRNAP ($factor) association with sigma factor\t1 $factors{RNAP}{$factor} + 1 $sigm --> 1 $rnap\tirreversible\tRecycling\n";
    }
    push @{$reactions{"Recycling"}}, "fmet_tRNA_cycle\tCharging of fmet-tRNA with methionine and formyl group\t1 fmet_tRNA + 1 $cpd_map{for} + 1 $cpd_map{Met} --> 1 fmet_tRNA_met\tirreversible\tRecycling\n";

    my $report = "Genetic Code is $assigned_code\n";
    my $reportObj = { "text_message"=>$report, "objects_created"=>[] };
    my $reportName = "build_me_model_report";

    my $time = time;
    my $metadata = $wsClient->save_objects({
	'workspace' => $workspace,
	'objects' => [{
	    type => 'KBaseReport.Report',
	    data => $reportObj,
	    name => $reportName,
	    'meta' => {},
	    'hidden' => 1,
	    'provenance' => $provenance
		      }]});

    $return = { 'report_name'=>$reportName, 'report_ref', $metadata->[0]->[6]."/".$metadata->[0]->[0]."/".$metadata->[0]->[4], 'reactions'=>\%reactions, 'compounds'=>\%compounds };

    #END build_me_model
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to build_me_model:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'build_me_model');
    }
    return($return);
}




=head2 version 

  $return = $obj->version()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
}

=head1 TYPES



=head2 ws_genome_id

=over 4



=item Description

The workspace ID for a ContigSet data object.
@id ws KBaseGenomes.Genome


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 workspace_name

=over 4



=item Description

A string representing a workspace name.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 MEModelingResult

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string


=end text

=back



=cut

1;
