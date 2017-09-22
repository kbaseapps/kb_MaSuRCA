# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import os
from pprint import pprint, pformat
from MaSuRCA.core.masurca_assembler import MaSuRCA_Assembler
#END_HEADER


class kb_MaSuRCA:
    '''
    Module Name:
    kb_MaSuRCA

    Module Description:
    Name of module: kb_MaSuRCA

This KBase module wraps the genome assembly software MaSuRCA(Maryland Super-Read Celera Assembler).
MaSuRCA 3.2.3


References:
https://academic.oup.com/bioinformatics/article/29/21/2669/195975/The-MaSuRCA-genome-assembler
https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btt476
ftp://ftp.genome.umd.edu/pub/MaSuRCA/latest/
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/kbaseapps/kb_MaSuRCA.git"
    GIT_COMMIT_HASH = "c67dbc6062d5cb4673d7a5058e5b6c7aecbfe032"

    #BEGIN_CLASS_HEADER
    # Class variables and functions can be defined in this block
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR

        # Any configuration parameters that are important should be parsed and
        # saved in the constructor.
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.work_folder = config['scratch']
        #END_CONSTRUCTOR
        pass


    def run_masurca(self, ctx, params):
        """
        Definition of run_masurca
        :param params: instance of type "masurcaParams" (Arguments for
           run_masurca *******for creating the sr_config.txt file******* 1.
           DATA consisting of 5 fields: 1)two_letter_prefix 2)mean 3)stdev
           4)fastq(.gz)_fwd_reads 5)fastq(.gz)_rev_reads. e.g., PE= pe 180 20
           /FULL_PATH/frag_1.fastq  /FULL_PATH/frag_2.fastq JUMP= sh 3600 200
           /FULL_PATH/short_1.fastq  /FULL_PATH/short_2.fastq #pacbio reads
           must be in a single fasta file! make sure you provide absolute
           path PACBIO=/FULL_PATH/pacbio.fa OTHER=/FULL_PATH/file.frg 2.
           PARAMETERS string GRAPH_KMER_SIZE - the k-mer size for deBruijn
           graph values between 25 and 127 are supported, 'auto' will compute
           the optimal size based on the read data and GC content bool
           USE_LINKING_MATES - set this to 1 for all Illumina-only
           assemblies; set this to 1 if you have less than 20x long reads
           (454, Sanger, Pacbio) and less than 50x CLONE coverage by
           Illumina, Sanger or 454 mate pairs; otherwise keep at 0 int
           LIMIT_JUMP_COVERAGE - this parameter is useful if you have too
           many Illumina jumping library mates. Typically set it to 60 for
           bacteria and 300 for the other organisms CA_PARAMETERS: these are
           the additional parameters to Celera Assembler.  do not worry about
           performance, number or processors or batch sizes -- these are
           computed automatically. float cgwErrorRate=0.15 - set
           cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for
           other organisms. int KMER_COUNT_THRESHOLD - minimum count k-mers
           used in error correction 1 means all k-mers are used.  one can
           increase to 2 if Illumina coverage >100 bool CLOSE_GAPS - whether
           to attempt to close gaps in scaffolds with Illumina data (1) or
           not (0) int NUM_THREADS - auto-detected number of cpus to use,
           mandatory int JF_SIZE  - this is mandatory jellyfish hash size --
           a safe value is estimated_genome_size*estimated_coverage (e.g.,
           2000000000) bool SOAP_ASSEMBLY - set this to 1 to use SOAPdenovo
           contigging/scaffolding module.  Assembly will be worse but will
           run faster. Useful for very large (>5Gbp) genomes string
           workspace_name - the name of the workspace from which to take
           input and store output. int hash_length - an odd integer (if even,
           it will be decremented) <= 31 string output_contigset_name - the
           name of the output contigset list<paired_end_lib> read_libraries -
           Illumina PairedEndLibrary files to assemble min_contig_length -
           integer to filter out contigs with length < min_contig_length from
           the MaSuRCA output. Default value is 500 (where 0 implies no
           filter). @optional min_contig_length @optional cov_cutoff
           @optional ins_length @optional read_trkg @optional amos_file
           @optional exp_cov @optional long_cov_cutoff) -> structure:
           parameter "workspace_name" of String, parameter "hash_length" of
           Long, parameter "read_libraries" of list of type "read_lib" (The
           workspace object name of a SingleEndLibrary or PairedEndLibrary
           file, whether of the KBaseAssembly or KBaseFile type.), parameter
           "output_contigset_name" of String, parameter "min_contig_length"
           of Long, parameter "cov_cutoff" of Double, parameter "ins_length"
           of Long, parameter "read_trkg" of type "bool" (A boolean - 0 for
           false, 1 for true. @range (0, 1)), parameter "amos_file" of type
           "bool" (A boolean - 0 for false, 1 for true. @range (0, 1)),
           parameter "exp_cov" of Double, parameter "long_cov_cutoff" of
           Double
        :returns: instance of type "masurcaResults" (Output parameter items
           for run_masurca report_name - the name of the KBaseReport.Report
           workspace object. report_ref - the workspace reference of the
           report.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_masurca
        self.log('Running run_masurca with params:\n' + pformat(params))
        for key, value in params.iteritems():
            if isinstance(value, basestring):
                params[key] = value.strip()

        masurca_assembler = MaSuRCA_Assembler(self.config, ctx.provenance())

        output = masurca_assembler.run_masurca(params)
        #END run_masurca

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_masurca return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
