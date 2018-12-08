# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import os
import json
import time
from MaSuRCA.core.masurca_assembler import MaSuRCA_Assembler
#END_HEADER


class kb_MaSuRCA:
    '''
    Module Name:
    kb_MaSuRCA

    Module Description:
    Name of module: MaSuRCA

This KBase module wraps the genome assembly software MaSuRCA(Maryland Super-Read Celera Assembler).
MaSuRCA 3.2.9


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
    VERSION = "1.1.0"
    GIT_URL = "https://github.com/kbaseapps/kb_MaSuRCA.git"
    GIT_COMMIT_HASH = "bb5563cc23532d6d78cf0d038c76e5d70319b13a"

    #BEGIN_CLASS_HEADER
    # Class variables and functions can be defined in this block

    def log(self, message, prefix_newline=False):
        print(('\n' if prefix_newline else '') +
              str(time.time()) + ': ' + str(message))
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR

        # Any configuration parameters that are important should be parsed and
        # saved in the constructor.
        self.config = config
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.config['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']
        #END_CONSTRUCTOR
        pass


    def run_masurca_assembler(self, ctx, params):
        """
        Definition of run_masurca_assembler
        :param params: instance of type "masurcaAssemblerParams" ->
           structure: parameter "workspace_name" of String, parameter
           "num_threads" of Long, parameter "jf_size" of Long, parameter
           "reads_libraries" of list of type "paired_readsParams" (parameter
           groups) -> structure: parameter "pe_id" of type "obj_ref" (An
           X/Y/Z style KBase object reference), parameter "pe_prefix" of
           String, parameter "pe_mean" of Long, parameter "pe_stdev" of Long,
           parameter "jump_libraries" of list of type "jump_readsParams" ->
           structure: parameter "jp_id" of type "obj_ref" (An X/Y/Z style
           KBase object reference), parameter "jp_prefix" of String,
           parameter "jp_mean" of Long, parameter "jp_stdev" of Long,
           parameter "pacbio_reads" of type "obj_ref" (An X/Y/Z style KBase
           object reference), parameter "nanopore_reads" of type "obj_ref"
           (An X/Y/Z style KBase object reference), parameter
           "other_frg_file" of String, parameter "graph_kmer_size" of String,
           parameter "use_linking_mates" of type "bool" (A boolean - 0 for
           false, 1 for true. @range (0, 1)), parameter "limit_jump_coverage"
           of Long, parameter "cgwErrorRate" of Double, parameter
           "kmer_count_threshold" of Long, parameter "close_gaps" of type
           "bool" (A boolean - 0 for false, 1 for true. @range (0, 1)),
           parameter "soap_assembly" of type "bool" (A boolean - 0 for false,
           1 for true. @range (0, 1)), parameter "do_homopolymer_trim" of
           type "bool" (A boolean - 0 for false, 1 for true. @range (0, 1)),
           parameter "output_contigset_name" of String, parameter
           "create_report" of type "bool" (A boolean - 0 for false, 1 for
           true. @range (0, 1))
        :returns: instance of type "masurcaResults" (Output parameter items
           for run_masurca and run_masurca_assembler report_name - the name
           of the KBaseReport.Report workspace object. report_ref - the
           workspace reference of the report.) -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_masurca_assembler
        self.log('Running run_masurca_assembler with params:\n{}'.format(
                 json.dumps(params, indent=1)))

        for key, value in params.iteritems():
            if isinstance(value, str):
                params[key] = value.strip()

        masurca_assembler = MaSuRCA_Assembler(self.config, ctx.provenance())

        output = masurca_assembler.run_masurca_assembler(params)
        #END run_masurca_assembler

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_masurca_assembler return value ' +
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
