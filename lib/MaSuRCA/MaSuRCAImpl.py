# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import os
from Bio import SeqIO
from pprint import pprint, pformat
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from KBaseReport.KBaseReportClient import KBaseReport
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
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

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

        #END_CONSTRUCTOR
        pass


    def run_masurca(self, ctx, params):
        """
        Definition of run_masurca
        :param params: instance of type "masurcaParams" (Arguments for
           run_masurca string workspace_name - the name of the workspace from
           which to take input and store output. int hash_length - an odd
           integer (if even, it will be decremented) <= 31 string
           output_contigset_name - the name of the output contigset
           list<paired_end_lib> read_libraries - Illumina PairedEndLibrary
           files to assemble min_contig_length - integer to filter out
           contigs with length < min_contig_length from the MaSuRCA output.
           Default value is 500 (where 0 implies no filter). @optional
           min_contig_length @optional cov_cutoff @optional ins_length
           @optional read_trkg @optional amos_file @optional exp_cov
           @optional long_cov_cutoff) -> structure: parameter
           "workspace_name" of String, parameter "hash_length" of Long,
           parameter "read_libraries" of list of type "read_lib" (The
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
