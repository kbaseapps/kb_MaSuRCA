import re
import time
import json
import os
import psutil
from pprint import pprint, pformat

from MaSuRCA.core.Program_Runner import Program_Runner
from DataFileUtil.DataFileUtilClient import DataFileUtil
from Workspace.WorkspaceClient import Workspace as Workspace
from KBaseReport.KBaseReportClient import KBaseReport
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from kb_quast.kb_quastClient import kb_quast


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))

class masurca_utils:
    MaSuRCA_VERSION = 'MaSuRCA-3.2.3'
    MaSuRCA_BIN = '/kb/module/' + MaSuRCA_VERSION + 'bin/masurca'
    PARAM_IN_WS = 'workspace_name'
    PARAM_IN_THREADN = 'num_threads'
    PARAM_IN_READS_LIBS = 'reads_libraries'
    PARAM_IN_JUMP_LIBS = 'jump_libraries'
    PARAM_IN_JF_SIZE = 'jf_size'
    PARAM_IN_CS_NAME = 'output_contigset_name'


    def __init__(self, scratch_dir, workspace_url, callback_url, provenance):
        self.workspace_url = workspace_url
        self.callback_url = callback_url
        self.au = AssemblyUtil(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url, service_ver='beta')
        self.scratch = scratch_dir
        self.working_dir = scratch_dir
        self.prog_runner = Program_Runner(self.MaSuRCA_BIN)
        self.provenance = provenance
        self.ws_client = Workspace(self.workspace_url)
        self.kbq = kb_quast(self.callbackURL)


    def validate_params(self, params):
        """
        validate_params: checks params passed to run_masurca_app method and set default values
        """
        log('Start validating run_masurca_app parameters')

        # check for mandatory parameters
        if params.get(self.PARAM_IN_WS, None) is None:
            raise ValueError(self.PARAM_IN_WS + ' parameter is required')
        if self.PARAM_IN_THREADN not in params:
            raise ValueError(self.PARAM_IN_THREADN + ' parameter is required')
        params[self.PARAM_IN_THREADN] = min(params.get(self.PARAM_IN_THREADN), psutil.cpu_count())

        if params.get(self.PARAM_IN_JF_SIZE, None) is None:
            raise ValueError(self.PARAM_IN_JF_SIZE + ' parameter is required')
        if self.PARAM_IN_READS_LIB not in params:
            raise ValueError(self.PARAM_IN_READS_LIB + ' parameter is required')
        if type(params[self.PARAM_IN_READS_LIB]) != list:
            raise ValueError(self.PARAM_IN_READS_LIB + ' must be a list')
        if (params.get(self.PARAM_IN_CS_NAME, None) is None or
                not valid_string(params[self.PARAM_IN_CS_NAME])):
            raise ValueError("Parameter output_contigset_name is required and must be a valid Workspace object string, "
                      "not {}".format(params.get(self.PARAM_IN_CS_NAME, None)))
        if ('pe_prefix' not in params):
            params['pe_prefix'] = 'pe'
        if ('pe_mean' not in params or type(params['pe_mean']) != int):
            params['pe_mean'] = 180
        if ('pe_stdv' not in params or type(params['pe_stdv']) != int):
            params['pe_stdv'] = 20

        return params

    def construct_masurca_config(self, params):
        # STEP 1: get the working folder housing the config.txt file and the masurca results
        out_folder = params['out_folder']
        out_dir = os.path.join(self.scratch, out_folder)
        self._mkdir_p(out_dir)

        wsname = params[self.PARAM_IN_WS]
        config_file_path = os.path.join(out_dir, 'config.txt')

        # STEP 2: retrieve the reads data from input parameter
        pe_reads_data = self._getKBReadsInfo(params[self.PARAM_IN_READS_LIBS])
        jp_reads_data = self._getKBReadsInfo(params[self.PARAM_IN_JUMP_LIBS])

        # STEP 3: construct and save the config.txt file for running masurca
        try:
            # STEP 3.1: replace the 'DATA...END' portion of the config_template.txt file 
            with open(config_file_path, 'w') as config_file:
                with open(os.path.join(os.path.dirname(__file__), 'config_template.txt'),
                          'r') as config_template_file:
                    config_template = config_template_file.read()
                    data_str = ''
                    if pe_reads_data:
                        data_str += 'PE= ' + params['pe_prefix'] + ' ' + str(params['pe_mean']) + ' ' + str(params['pe_stdv']) + ' ' + pe_reads_data['fwd_file'] + ' ' + pe_reads_data['reverse_file']
                    if jp_reads_data:
                        if ('jp_mean' not in params or type(params['jp_mean']) != int):
                            params['jp_mean'] = 3600
                        if ('pe_stdv' not in params or type(params['jp_stdv']) != int):
                            params['pe_stdv'] = 200
                        if data_str != '':
                            data_str += '\n'
                        data_str += 'JUMP= ' + params['jp_prefix'] + ' ' + str(params['jp_mean']) + ' ' + str(params['jp_stdv']) + ' ' + jp_reads_data['fwd_file'] + ' ' + jp_reads_data['reverse_file']

                    begin_patn1 = "DATA\n"
                    end_patn1 = "END\nPARAMETERS\n"
                    config_template = self._replaceSectionText(config_template, begin_patn1, end_patn1, data_str)
                    config_file.write(config_template)

            # STEP 3.2: replace the 'PARAMETERS...END' portion of the config_file file saved in STEP 3.1
            with open(config_file_path, 'r') as previous_config_file:
                previous_config = previous_config_file.read()
                param_str = ''
                if params['graph_kmer_size']:
                    param_str += 'GRAPH_KMER_SIZE=' + str(params['graph_kmer_size'])
                if params['use_linking_mates'] == 1:
                    param_str += '\nUSE_LINKING_MATES=1'
                else:
                    param_str += '\nUSE_LINKING_MATES=0'
                if params['limit_jump_coverage']:
                    param_str += '\nLIMIT_JUMP_COVERAGE = ' + str(params['limit_jump_coverage'])
                if params['cgwErrorRate']:
                    param_str += '\nCA_PARAMETERS = cgwErrorRate=' + str(params['cgwErrorRate'])
                if params['num_threads']:
                    param_str += '\nNUM_THREADS=' + str(params['num_threads'])
                if params['jf_size']:
                    param_str += '\nJF_SIZE=' + str(params['jf_size'])
                if params['do_homopolymer_trim'] == 1:
                    param_str += '\nDO_HOMOPOLYMER_TRIM==1'
                else:
                    param_str += '\nDO_HOMOPOLYMER_TRIM==0'

                begin_patn2 = "PARAMETERS\n"
                end_patn2 = "END\n"
                param_str = begin_patn + param_str + + '\n' + end_patn
                final_config = self._replaceSectionText(previous_config, begin_patn2, end_patn2, param_str)

            with open(config_file_path, 'w') as config_file:
                config_file.write(final_config)
        except IOError as ioerr:
            log('Creation of the config.txt file raised error:\n')
            pprint(ioerr)
            return ''
        else:
            return config_file_path

    def generate_assemble_script(self, config_file):
        exit_code = 1
        if config_file != '':
            f_dir, f_nm = os.path.split(config_file)
            m_cmd = [self.MaSuRCA_BIN]
            m_cmd.append(config_file)
            exit_code = self.prog_runner.run(m_cmd, f_dir)

        if exit_code == 0:
            return os.path.join(f_dir, 'assemble.sh')
        else:
            return ''

    def run_assemble(self, asmbl_file):
        exit_code = 1
        if asmbl_file != '':
            f_dir, f_nm = os.path.split(asmbl_file)
            a_cmd = []
            a_cmd.append('./' + asmbl_file)
            exit_code = self.prog_runner.run(a_cmd, f_dir)

        return exit_code

    def _replaceSectionText(self, orig_txt, begin_patn, end_patn, repl_txt):
        """
        replace a section of text of orig_txt between lines begin-patn and end-patn with repl_text
        examples of parameters:
            begin_patn1 = "DATA\n"
            begin_patn2 = "PARAMETERS\n"
            end_patn1 = "END\nPARAMETERS\n"
            end_patn2 = "END\n"
            repl_txt1 = 'PE= pe 180 20 /kb/module/work/testReads/small.forward.fq /kb/module/work/testReads/small.reverse.fq\n'
            repl_txt2 = 'GRAPH_KMER_SIZE=auto\nUSE_LINKING_MATES=1\nLIMIT_JUMP_COVERAGE = 60\nCA_PARAMETERS = cgwErrorRate=0.15\nNUM_THREADS= 64\nJF_SIZE=100000000\nDO_HOMOPOLYMER_TRIM=0\n' 
        """
        if repl_txt != '':
            # create regular expression pattern
            repl = re.compile(begin_patn + '.*?' + end_patn, re.DOTALL)

            repl_txt = begin_patn + repl_txt + '\n' + end_patn
            # replace the text between begin_patn and end_patn with repl_txt
            txt_replaced = repl.sub(repl_txt, orig_txt)
            #pprint(txt_replaced)
            return txt_replaced
        else:
            return orig_txt


    def _getKBReadsInfo(self, wsname, reads_refs):
        """
        _getKBReadsInfo--from a set of given KBase reads refs, fetch the corresponding reads info
        and return the results in the following structure:
        reads_date = {
                'fwd_file': path_to_fastq_file,
                'type': reads_type, #('interleaved', 'paired', or 'single'
                'rev_file': path_to_fastq_file, #only if paired end
                'seq_tech': sequencing_tech,
                'reads_ref': KBase object ref for downstream convenience
        }
        """
        obj_ids = []
        for r in reads_refs:
            obj_ids.append({'ref': r if '/' in r else (wsname + '/' + r)})

        ws = workspaceService(self.workspaceURL, token=token)
        ws_info = ws.get_object_info_new({'objects': obj_ids})
        reads_params = []

        reftoname = {}
        for wsi, oid in zip(ws_info, obj_ids):
            ref = oid['ref']
            reads_params.append(ref)
            obj_name = wsi[1]
            reftoname[ref] = wsi[7] + '/' + obj_name

        typeerr = ('Supported types: KBaseFile.SingleEndLibrary ' +
                   'KBaseFile.PairedEndLibrary ' +
                   'KBaseAssembly.SingleEndLibrary ' +
                   'KBaseAssembly.PairedEndLibrary')
        try:
            readcli = ReadsUtils(self.callbackURL, token=token)
            reads = readcli.download_reads({'read_libraries': reads_params})['files']
        except ServerError as se:
            self.log('logging stacktrace from dynamic client error')
            self.log(se.data)
            if typeerr in se.message:
                prefix = se.message.split('.')[0]
                raise ValueError(
                    prefix + '. Only the types ' +
                    'KBaseAssembly.SingleEndLibrary ' +
                    'KBaseAssembly.PairedEndLibrary ' +
                    'KBaseFile.SingleEndLibrary ' +
                    'and KBaseFile.PairedEndLibrary are supported')
            else:
                raise

        self.log('Got reads data from converter:\n' + pformat(reads))

        reads_data = []
        for ref in reads:
            reads_name = reftoname[ref]
            f = reads[ref]['files']
            seq_tech = reads[ref]["sequencing_tech"]
            if f['type'] == 'interleaved':
                reads_data.append({'fwd_file': f['fwd'], 'type':'interleaved',
                                   'seq_tech': seq_tech, 'reads_ref': ref})
            elif f['type'] == 'paired':
                reads_data.append({'fwd_file': f['fwd'], 'rev_file': f['rev'],
                                   'type':'paired', 'seq_tech': seq_tech, 'reads_ref': ref})
            elif f['type'] == 'single':
                reads_data.append({'fwd_file': f['fwd'], 'type':'single',
                                   'seq_tech': seq_tech, 'reads_ref': ref})
            else:
                raise ValueError('Something is very wrong with read lib' + reads_name)

        return reads_data

    def _mkdir_p(self, dir):
        """
        _mkdir_p: make directory for given path
        """
        log('Creating a new dir: ' + dir)
        if not dir:
            return
        if not os.path.exists(dir):
            os.makedirs(dir)
        else:
            log('{} has existed, so skip creating.'.format(dir))


