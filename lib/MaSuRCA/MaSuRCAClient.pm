package MaSuRCA::MaSuRCAClient;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

MaSuRCA::MaSuRCAClient

=head1 DESCRIPTION


Name of module: MaSuRCA

This KBase module wraps the genome assembly software MaSuRCA(Maryland Super-Read Celera Assembler).
MaSuRCA 3.2.9


References:
https://academic.oup.com/bioinformatics/article/29/21/2669/195975/The-MaSuRCA-genome-assembler
https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btt476
ftp://ftp.genome.umd.edu/pub/MaSuRCA/latest/


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => MaSuRCA::MaSuRCAClient::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my %arg_hash2 = @args;
	if (exists $arg_hash2{"token"}) {
	    $self->{token} = $arg_hash2{"token"};
	} elsif (exists $arg_hash2{"user_id"}) {
	    my $token = Bio::KBase::AuthToken->new(@args);
	    if (!$token->error_message) {
	        $self->{token} = $token->token;
	    }
	}
	
	if (exists $self->{token})
	{
	    $self->{client}->{token} = $self->{token};
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 run_masurca_assembler

  $output = $obj->run_masurca_assembler($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_MaSuRCA.masurcaAssemblerParams
$output is a kb_MaSuRCA.masurcaResults
masurcaAssemblerParams is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a string
	num_threads has a value which is an int
	jf_size has a value which is an int
	reads_libraries has a value which is a reference to a list where each element is a kb_MaSuRCA.paired_readsParams
	jump_libraries has a value which is a reference to a list where each element is a kb_MaSuRCA.jump_readsParams
	pacbio_reads has a value which is a kb_MaSuRCA.obj_ref
	nanopore_reads has a value which is a kb_MaSuRCA.obj_ref
	other_frg_file has a value which is a string
	graph_kmer_size has a value which is a string
	use_linking_mates has a value which is a kb_MaSuRCA.bool
	dna_source has a value which is a string
	kmer_count_threshold has a value which is an int
	close_gaps has a value which is a kb_MaSuRCA.bool
	soap_assembly has a value which is a kb_MaSuRCA.bool
	do_homopolymer_trim has a value which is a kb_MaSuRCA.bool
	output_contigset_name has a value which is a string
	create_report has a value which is a kb_MaSuRCA.bool
paired_readsParams is a reference to a hash where the following keys are defined:
	pe_id has a value which is a kb_MaSuRCA.obj_ref
	pe_prefix has a value which is a string
	pe_mean has a value which is an int
	pe_stdev has a value which is an int
obj_ref is a string
jump_readsParams is a reference to a hash where the following keys are defined:
	jp_id has a value which is a kb_MaSuRCA.obj_ref
	jp_prefix has a value which is a string
	jp_mean has a value which is an int
	jp_stdev has a value which is an int
bool is an int
masurcaResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string

</pre>

=end html

=begin text

$params is a kb_MaSuRCA.masurcaAssemblerParams
$output is a kb_MaSuRCA.masurcaResults
masurcaAssemblerParams is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a string
	num_threads has a value which is an int
	jf_size has a value which is an int
	reads_libraries has a value which is a reference to a list where each element is a kb_MaSuRCA.paired_readsParams
	jump_libraries has a value which is a reference to a list where each element is a kb_MaSuRCA.jump_readsParams
	pacbio_reads has a value which is a kb_MaSuRCA.obj_ref
	nanopore_reads has a value which is a kb_MaSuRCA.obj_ref
	other_frg_file has a value which is a string
	graph_kmer_size has a value which is a string
	use_linking_mates has a value which is a kb_MaSuRCA.bool
	dna_source has a value which is a string
	kmer_count_threshold has a value which is an int
	close_gaps has a value which is a kb_MaSuRCA.bool
	soap_assembly has a value which is a kb_MaSuRCA.bool
	do_homopolymer_trim has a value which is a kb_MaSuRCA.bool
	output_contigset_name has a value which is a string
	create_report has a value which is a kb_MaSuRCA.bool
paired_readsParams is a reference to a hash where the following keys are defined:
	pe_id has a value which is a kb_MaSuRCA.obj_ref
	pe_prefix has a value which is a string
	pe_mean has a value which is an int
	pe_stdev has a value which is an int
obj_ref is a string
jump_readsParams is a reference to a hash where the following keys are defined:
	jp_id has a value which is a kb_MaSuRCA.obj_ref
	jp_prefix has a value which is a string
	jp_mean has a value which is an int
	jp_stdev has a value which is an int
bool is an int
masurcaResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string


=end text

=item Description

Definition of run_masurca_assembler

=back

=cut

 sub run_masurca_assembler
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function run_masurca_assembler (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to run_masurca_assembler:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'run_masurca_assembler');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_MaSuRCA.run_masurca_assembler",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'run_masurca_assembler',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method run_masurca_assembler",
					    status_line => $self->{client}->status_line,
					    method_name => 'run_masurca_assembler',
				       );
    }
}
 
  
sub status
{
    my($self, @args) = @_;
    if ((my $n = @args) != 0) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function status (received $n, expecting 0)");
    }
    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
        method => "kb_MaSuRCA.status",
        params => \@args,
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => 'status',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
                          );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method status",
                        status_line => $self->{client}->status_line,
                        method_name => 'status',
                       );
    }
}
   

sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_MaSuRCA.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'run_masurca_assembler',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method run_masurca_assembler",
            status_line => $self->{client}->status_line,
            method_name => 'run_masurca_assembler',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for MaSuRCA::MaSuRCAClient\n";
    }
    if ($sMajor == 0) {
        warn "MaSuRCA::MaSuRCAClient version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 assembly_ref

=over 4



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



=head2 bool

=over 4



=item Description

A boolean - 0 for false, 1 for true.
@range (0, 1)


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 obj_ref

=over 4



=item Description

An X/Y/Z style KBase object reference


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



=head2 paired_readsParams

=over 4



=item Description

parameter groups


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
pe_id has a value which is a kb_MaSuRCA.obj_ref
pe_prefix has a value which is a string
pe_mean has a value which is an int
pe_stdev has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
pe_id has a value which is a kb_MaSuRCA.obj_ref
pe_prefix has a value which is a string
pe_mean has a value which is an int
pe_stdev has a value which is an int


=end text

=back



=head2 jump_readsParams

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
jp_id has a value which is a kb_MaSuRCA.obj_ref
jp_prefix has a value which is a string
jp_mean has a value which is an int
jp_stdev has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
jp_id has a value which is a kb_MaSuRCA.obj_ref
jp_prefix has a value which is a string
jp_mean has a value which is an int
jp_stdev has a value which is an int


=end text

=back



=head2 masurcaAssemblerParams

=over 4



=item Description

Arguments for run_masurca

*******for creating the sr_config.txt file*******
1. DATA
consisting of 5 fields: 1)two_letter_prefix 2)mean 3)stdev 4)fastq(.gz)_fwd_reads 5)fastq(.gz)_rev_reads.
e.g.,
        PE= pe 180 20  /FULL_PATH/frag_1.fastq  /FULL_PATH/frag_2.fastq
        JUMP= sh 3600 200  /FULL_PATH/short_1.fastq  /FULL_PATH/short_2.fastq
        #pacbio OR nanopore reads must be in a single fasta or fastq file with absolute path, can be gzipped
        #if you have both types of reads supply them both as NANOPORE type
        PACBIO=/FULL_PATH/pacbio.fa
        NANOPORE=/FULL_PATH/nanopore.fa
        OTHER=/FULL_PATH/file.frg

2. PARAMETERS
string graph_kmer_size - the k-mer size for deBruijn graph values between 25 and 127 are supported, 'auto' will compute the optimal size based on the read data and GC content
bool use_linking_mates - set this to 1 for all Illumina-only assemblies; set this to 1 if you have less than 20x long reads (454, Sanger, Pacbio) and less than 50x CLONE coverage by Illumina, Sanger or 454 mate pairs; otherwise keep at 0
string dna_source - indicate 'bacteria' or 'other organisms' for setting limit_jump_coverage and cgwErrorRate values
int limit_jump_coverage - this parameter is useful if you have too many Illumina jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms
CA_PARAMETERS: these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. 
float cgwErrorRate=0.15 - set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
int kmer_count_threshold - minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if Illumina coverage >100
bool close_gaps - whether to attempt to close gaps in scaffolds with Illumina data (1) or not (0)
int num_threads - auto-detected number of cpus to use, mandatory
int jf_size  - this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*estimated_coverage (e.g., 2000000000)
bool SOAP_ASSEMBLY - set this to 1 to use SOAPdenovo contigging/scaffolding module.  Assembly will be worse but will run faster. Useful for very large (>5Gbp) genomes
bool do_homopolymer_trim - specifies if we do (1) or do not (0) want to trim long runs of homopolymers 

string workspace_name - the name of the workspace from which to take input and store output.
string output_contigset_name - the name of the output contigset
list<paired_readsParams> read_libraries - Illumina PairedEndLibrary files to assemble

@optional jump_libraries
@optional pacbio_reads
@optional other_frg_file
@optional graph_kmer_size
@optional use_linking_mates
@optional dna_source
@optional kmer_count_threshold
@optional close_gaps
@optional soap_assembly
@optional do_homopolymer_trim


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a string
num_threads has a value which is an int
jf_size has a value which is an int
reads_libraries has a value which is a reference to a list where each element is a kb_MaSuRCA.paired_readsParams
jump_libraries has a value which is a reference to a list where each element is a kb_MaSuRCA.jump_readsParams
pacbio_reads has a value which is a kb_MaSuRCA.obj_ref
nanopore_reads has a value which is a kb_MaSuRCA.obj_ref
other_frg_file has a value which is a string
graph_kmer_size has a value which is a string
use_linking_mates has a value which is a kb_MaSuRCA.bool
dna_source has a value which is a string
kmer_count_threshold has a value which is an int
close_gaps has a value which is a kb_MaSuRCA.bool
soap_assembly has a value which is a kb_MaSuRCA.bool
do_homopolymer_trim has a value which is a kb_MaSuRCA.bool
output_contigset_name has a value which is a string
create_report has a value which is a kb_MaSuRCA.bool

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a string
num_threads has a value which is an int
jf_size has a value which is an int
reads_libraries has a value which is a reference to a list where each element is a kb_MaSuRCA.paired_readsParams
jump_libraries has a value which is a reference to a list where each element is a kb_MaSuRCA.jump_readsParams
pacbio_reads has a value which is a kb_MaSuRCA.obj_ref
nanopore_reads has a value which is a kb_MaSuRCA.obj_ref
other_frg_file has a value which is a string
graph_kmer_size has a value which is a string
use_linking_mates has a value which is a kb_MaSuRCA.bool
dna_source has a value which is a string
kmer_count_threshold has a value which is an int
close_gaps has a value which is a kb_MaSuRCA.bool
soap_assembly has a value which is a kb_MaSuRCA.bool
do_homopolymer_trim has a value which is a kb_MaSuRCA.bool
output_contigset_name has a value which is a string
create_report has a value which is a kb_MaSuRCA.bool


=end text

=back



=head2 masurcaResults

=over 4



=item Description

Output parameter items for run_masurca_assembler
report_name - the name of the KBaseReport.Report workspace object.
report_ref - the workspace reference of the report.


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

package MaSuRCA::MaSuRCAClient::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
