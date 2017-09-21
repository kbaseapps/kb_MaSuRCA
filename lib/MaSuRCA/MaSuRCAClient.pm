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


Name of module: kb_MaSuRCA

This KBase module wraps the genome assembly software MaSuRCA(Maryland Super-Read Celera Assembler).
MaSuRCA 3.2.3


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




=head2 run_masurca

  $output = $obj->run_masurca($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_MaSuRCA.masurcaParams
$output is a kb_MaSuRCA.masurcaResults
masurcaParams is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a string
	hash_length has a value which is an int
	read_libraries has a value which is a reference to a list where each element is a kb_MaSuRCA.read_lib
	output_contigset_name has a value which is a string
	min_contig_length has a value which is an int
	cov_cutoff has a value which is a float
	ins_length has a value which is an int
	read_trkg has a value which is a kb_MaSuRCA.bool
	amos_file has a value which is a kb_MaSuRCA.bool
	exp_cov has a value which is a float
	long_cov_cutoff has a value which is a float
read_lib is a string
bool is an int
masurcaResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string

</pre>

=end html

=begin text

$params is a kb_MaSuRCA.masurcaParams
$output is a kb_MaSuRCA.masurcaResults
masurcaParams is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a string
	hash_length has a value which is an int
	read_libraries has a value which is a reference to a list where each element is a kb_MaSuRCA.read_lib
	output_contigset_name has a value which is a string
	min_contig_length has a value which is an int
	cov_cutoff has a value which is a float
	ins_length has a value which is an int
	read_trkg has a value which is a kb_MaSuRCA.bool
	amos_file has a value which is a kb_MaSuRCA.bool
	exp_cov has a value which is a float
	long_cov_cutoff has a value which is a float
read_lib is a string
bool is an int
masurcaResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string


=end text

=item Description

Definition of run_masurca

=back

=cut

 sub run_masurca
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function run_masurca (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to run_masurca:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'run_masurca');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_MaSuRCA.run_masurca",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'run_masurca',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method run_masurca",
					    status_line => $self->{client}->status_line,
					    method_name => 'run_masurca',
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
                method_name => 'run_masurca',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method run_masurca",
            status_line => $self->{client}->status_line,
            method_name => 'run_masurca',
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



=head2 read_lib

=over 4



=item Description

The workspace object name of a SingleEndLibrary or PairedEndLibrary file, whether of the
KBaseAssembly or KBaseFile type.


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



=head2 masurcaParams

=over 4



=item Description

Arguments for run_masurca

string workspace_name - the name of the workspace from which to take input and store output.
int hash_length - an odd integer (if even, it will be decremented) <= 31
string output_contigset_name - the name of the output contigset
list<paired_end_lib> read_libraries - Illumina PairedEndLibrary files to assemble
min_contig_length - integer to filter out contigs with length < min_contig_length
             from the MaSuRCA output. Default value is 500 (where 0 implies no filter).

@optional min_contig_length
@optional cov_cutoff
@optional ins_length
@optional read_trkg
@optional amos_file
@optional exp_cov
@optional long_cov_cutoff


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a string
hash_length has a value which is an int
read_libraries has a value which is a reference to a list where each element is a kb_MaSuRCA.read_lib
output_contigset_name has a value which is a string
min_contig_length has a value which is an int
cov_cutoff has a value which is a float
ins_length has a value which is an int
read_trkg has a value which is a kb_MaSuRCA.bool
amos_file has a value which is a kb_MaSuRCA.bool
exp_cov has a value which is a float
long_cov_cutoff has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a string
hash_length has a value which is an int
read_libraries has a value which is a reference to a list where each element is a kb_MaSuRCA.read_lib
output_contigset_name has a value which is a string
min_contig_length has a value which is an int
cov_cutoff has a value which is a float
ins_length has a value which is an int
read_trkg has a value which is a kb_MaSuRCA.bool
amos_file has a value which is a kb_MaSuRCA.bool
exp_cov has a value which is a float
long_cov_cutoff has a value which is a float


=end text

=back



=head2 masurcaResults

=over 4



=item Description

Output parameter items for run_masurca

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
