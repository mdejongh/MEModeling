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
	$genome = $wsClient->get_objects([{ref=>$genome_id}])->[0]{data};
	push @{$provenance->[0]->{'input_ws_objects'}}, $genome_id;
    };
    if ($@) {
	die "Error loading genome:\n".$@;
    }

    my $assigned_code = $genome->{genetic_code};

    my $report = "Genetic Code is $assigned_code\n";
    my $reportObj = { "text_message"=>$report };
    my $reportName = "build_me_model_report";

    my $time = time;
    my $metadata = $wsClient->save_objects({
	'id' => $reportName.$time,
	'objects' => [{
	    type => 'KBaseReport.Report',
	    data => $reportObj,
	    name => $reportName,
	    'meta' => {},
	    'hidden' => 1,
	    'provenance' => $provenance
		      }]});

    $return = { 'report_name'=>$reportName, 'report_ref', $metadata->[0]->[6]."/".$metadata->[0]->[0]."/".$metadata->[0]->[4] };

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
