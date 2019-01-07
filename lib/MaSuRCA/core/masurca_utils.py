# -*- coding: utf-8 -*-
import re
import time
import json
import os
inport errno
import numpy as np
import zipfile
import subprocess
from pprint import pprint, pformat
import codecs
import uuid
import copy

from MaSuRCA.core.Program_Runner import Program_Runner
from Workspace.WorkspaceClient import Workspace as Workspace
from KBaseReport.KBaseReportClient import KBaseReport
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from kb_quast.kb_quastClient import kb_quast
from kb_quast.baseclient import ServerError as QUASTError
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from ReadsUtils.baseclient import ServerError


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


def _mkdir_p(path):
    """
    _mkdir_p: make directory for given path
    """
    if not path:
        return
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


class masurca_utils:
    MaSuRCA_VERSION = 'MaSuRCA-3.2.9'
    MaSuRCA_BIN = '/kb/module/' + MaSuRCA_VERSION + '/bin/masurca'
    PARAM_IN_WS = 'workspace_name'
    PARAM_IN_THREADN = 'num_threads'
    PARAM_IN_READS_LIBS = 'reads_libraries'
    PARAM_IN_JUMP_LIBS = 'jump_libraries'
    PARAM_IN_JF_SIZE = 'jf_size'
    PARAM_IN_CS_NAME = 'output_contigset_name'

    def __init__(self, prj_dir, config):
        self.workspace_url = config['workspace-url']
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        if 'shock-url' in config:
            self.shock_url = config['shock-url']
        if 'handle-service-url' in config:
            self.handle_url = config['handle-service-url']

        self.ws_client = Workspace(self.workspace_url, token=self.token)
        self.ru = ReadsUtils(self.callback_url, token=self.token)
        self.au = AssemblyUtil(self.callback_url, token=self.token)
        self.kbr = KBaseReport(self.callback_url)
        self.kbq = kb_quast(self.callback_url)
        self.proj_dir = prj_dir
        self.prog_runner = Program_Runner(self.MaSuRCA_BIN, self.proj_dir)

    def _has_long_reads(self, params):
        """
        _has_long_reads: check if a long reads input exists in the parameters
        """
        return (params.get('pacbio_reads', None) or
                params.get('nanopore_reads', None) or
                params.get('other_frg_file', None))

    def _get_data_portion(self, pe_reads_data, jp_reads_data=None, pacbio_reads_file='',
                          nanopore_reads_file='', other_frg_file=''):
        """
        _get_data_portion: build the 'DATA...END' portion for the config.txt file
        """
        data_str = ''
        if pe_reads_data:
            # log('PE reads data details:\n{}'.format(json.dumps(pe_reads_data, indent=1)))
            for pe in pe_reads_data:
                if data_str != '':
                    data_str += '\n'
                data_str += 'PE= ' + pe['pe_prefix'] + ' ' + str(pe['pe_mean']) + ' ' + \
                            str(pe['pe_stdev']) + ' ' + pe['fwd_file']
                if pe.get('rev_file', None):
                    data_str += ' ' + pe['rev_file']

        if jp_reads_data:
            # log('JUMP reads data details:\n{}'.format(json.dumps(jp_reads_data, indent=1)))
            for jp in jp_reads_data:
                if data_str != '':
                    data_str += '\n'
                data_str += 'JUMP= ' + jp['jp_prefix'] + ' ' + str(jp['jp_mean']) + ' ' + \
                            str(jp['jp_stdev']) + ' ' + jp['fwd_file']
                if jp.get('rev_file', None):
                    data_str += ' ' + jp['rev_file']

        # Adding the pacbio_reads
        # Note that pcbio reads must be in a single fasta file!
        # For example:
        # data_str +='\nPACBIO= /pool/genomics/frandsenp/masurca/PacBio/pacbio_reads.fasta'
        # ***if you have both types of reads supply them both as NANOPORE type***
        if pacbio_reads_file != '':
            if data_str != '':
                data_str += '\n'
            if nanopore_reads_file != '':
                data_str += 'NANOPORE=' + pacbio_reads_file
            else:
                data_str += 'PACBIO=' + pacbio_reads_file

        # Adding the nanopore_reads and note that nanopore reads must be in a single fasta file!
        # For example:
        # data_str +='\nNANOPORE= /pool/genomics/frandsenp/masurca/NanoPore/nanopore_reads.fasta'
        if nanopore_reads_file != '':
            if data_str != '':
                data_str += '\n'
            data_str += 'NANOPORE= ' + nanopore_reads_file

        # Adding the other_frg_file inputs if any
        # any OTHER sequence data (454, Sanger, Ion torrent, etc) must be first converted into 
        # Celera Assembler compatible .frg file
        # (see http://wgsassembler.sourceforge.com) and supplied as OTHER=file.frg
        if other_frg_file != '':
            if data_str != '':
                data_str += '\n'
            data_str += 'OTHER=' + other_frg_file

        return data_str

    def _get_parameters_portion(self, params):
        """
        build the 'PARAMETERS...END' portion for the config.txt file
        """
        # set the default parameters as suggested in the example configuration file
        param_str = ("EXTEND_JUMP_READS=0\nUSE_GRID=0\nGRID_QUEUE=all.q\nGRID_BATCH_SIZE" +
                     "=300000000\nLHE_COVERAGE=25\nMEGA_READS_ONE_PASS=0")
        if (params.get('graph_kmer_size', None) and
                type(params['graph_kmer_size']) == int):
            if param_str != '':
                param_str += '\n'
            param_str += 'GRAPH_KMER_SIZE=' + str(params['graph_kmer_size'])
        else:
            if param_str != '':
                param_str += '\n'
            param_str += 'GRAPH_KMER_SIZE=auto'
        if params.get('use_linking_mates', None):
            if param_str != '':
                param_str += '\n'
            if params['use_linking_mates'] == 1 and not self._has_long_reads(params):
                param_str += 'USE_LINKING_MATES=1'
            else:
                param_str += 'USE_LINKING_MATES=0'
        if params.get('limit_jump_coverage', None):
            if param_str != '':
                param_str += '\n'
            param_str += 'LIMIT_JUMP_COVERAGE = ' + str(params['limit_jump_coverage'])
        if params.get('cgwErrorRate', None):
            if param_str != '':
                param_str += '\n'
            param_str += 'CA_PARAMETERS = cgwErrorRate=' + str(params['cgwErrorRate'])
        if params.get(self.PARAM_IN_THREADN, None):
            if param_str != '':
                param_str += '\n'
            param_str += 'NUM_THREADS = ' + str(params[self.PARAM_IN_THREADN])
        if params.get('jf_size', None):
            if param_str != '':
                param_str += '\n'
            param_str += 'JF_SIZE=' + str(params['jf_size'])
        if params.get('kmer_count_threshold', None):
            if param_str != '':
                param_str += '\n'
            param_str += 'KMER_COUNT_THRESHOLD=' + str(params['kmer_count_threshold'])
        if params.get('do_homopolymer_trim', None):
            if param_str != '':
                param_str += '\n'
            if params['do_homopolymer_trim'] == 1:
                param_str += 'DO_HOMOPOLYMER_TRIM=1'
            else:
                param_str += 'DO_HOMOPOLYMER_TRIM=0'
        if params.get('close_gaps', None):
            if param_str != '':
                param_str += '\n'
            if params['close_gaps'] == 1:
                param_str += 'CLOSE_GAPS=1'
            else:
                param_str += 'CLOSE_GAPS=0'
        if params.get('soap_assembly', None):
            if param_str != '':
                param_str += '\n'
            if params['soap_assembly'] == 1:
                param_str += 'SOAP_ASSEMBLY=1'
            else:
                param_str += 'SOAP_ASSEMBLY=0'
        return param_str

    def _replaceSectionText(self, orig_txt, begin_patn, end_patn, repl_txt):
        """
        replace a section of text of orig_txt between lines begin-patn and end-patn with repl_text
        examples of parameters:
            begin_patn1 = "DATA\n"
            begin_patn2 = "PARAMETERS\n"
            end_patn1 = "END\nPARAMETERS\n"
            end_patn2 = "END\n"
            repl_txt1 = 'PE= pe 500 50 /kb/module/work/testReads/small.forward.fq /kb/module/work/testReads/small.reverse.fq\n'
            repl_txt2 = 'GRAPH_KMER_SIZE=auto\nUSE_LINKING_MATES=1\nLIMIT_JUMP_COVERAGE = 60\nCA_PARAMETERS = cgwErrorRate=0.15\nNUM_THREADS= 64\nJF_SIZE=100000000\nDO_HOMOPOLYMER_TRIM=0\n'
        """
        if repl_txt != '':
            # create regular expression pattern
            repl = re.compile(begin_patn + '.*?' + end_patn, re.DOTALL)
            repl_txt = begin_patn + repl_txt + '\n' + end_patn
            # replace the text between begin_patn and end_patn with repl_txt
            txt_replaced = repl.sub(repl_txt, orig_txt)
            # pprint(txt_replaced)
            return txt_replaced
        else:
            return orig_txt

    def _unique_prefix_check(self, pfix, refs):
        prefix_lookup = {}
        for ref in refs:
            pre = ref[pfix][0:2]
            if pre not in prefix_lookup:
                prefix_lookup[pre] = 1
            else:
                raise ValueError('The first two letters in \'' + ref + '\' has been used.')

    def _getReadsInfo_PE(self, input_params):
        """
        _getReadsInfo_PE--from a list of paired_readsParams structures fetches the corresponding
        reads info with the paired_readsParams[pe_id]
        returns a list of reads data in the following structure:
        reads_data = {
                'fwd_file': path_to_fastq_file,
                'pe_prefix': the two-letter prefix for the reads library,
                'pe_mean': the average reads length for the reads library,
                'pe_stdev': the standard deviation for the reads library,
                'type': reads_type, #('interleaved', 'paired', or 'single'
                'seq_tech': sequencing_tech,
                'reads_ref': KBase object ref for downstream convenience,
                'reads_name': KBase object name for downstream convenience,
                'rev_file': path_to_fastq_file, #only if paired end
        }
        """
        rds_params = copy.deepcopy(input_params)
        wsname = rds_params[self.PARAM_IN_WS]
        rds_refs = []
        rds_data = []
        self._unique_prefix_check('pe_prefix', rds_params[self.PARAM_IN_READS_LIBS])

        # reads_libraries grouped params
        if 'reads_libraries' in rds_params and rds_params['reads_libraries']:
            for rds_lib in rds_params['reads_libraries']:
                if 'pe_id' in rds_lib:
                    rds_refs.append(rds_lib['pe_id'])
            rds_data = self._getKBReadsInfo(wsname, rds_refs)

            for rds_lib in rds_params['reads_libraries']:
                i = 0
                for rds in rds_data:
                    i += 1
                    if ('pe_id' in rds_lib and rds_lib['pe_id'] == rds['reads_ref']):
                        if 'pe_prefix' in rds_lib:
                            rds['pe_prefix'] = rds_lib['pe_prefix'][:1] + str(i)
                        else:
                            raise ValueError("Parameter pe_prefix is required for reads {}".format(
                                rds[reads_ref]))
                        if 'pe_mean' in rds_lib:
                            rds['pe_mean'] = rds_lib['pe_mean']
                        else:
                            raise ValueError("Parameter pe_mean is required for reads {}".format(
                                rds[reads_ref]))
                        if 'pe_stdev' in rds_lib:
                            rds['pe_stdev'] = rds_lib['pe_stdev']
                        else:
                            raise ValueError("Parameter pe_stdev is required for reads {}".format(
                                rds[reads_ref]))
        else:
            raise ValueError("Parameter {} is required for reads {}".format('reads_libraries'))
        return rds_data

    def _getReadsInfo_JP(self, input_params):
        """
        _getReadsInfo_JP--from a list of jump_readsParams structures fetches the corresponding
        reads info with the paired_readsParams[pe_id]
        returns a list of reads data in the following structure:
        reads_data = {
                'fwd_file': path_to_fastq_file,
                'jp_prefix': the two-letter prefix for the reads library,
                'jp_mean': the average reads length for the reads library,
                'jp_stdev': the standard deviation for the reads library,
                'type': reads_type, #('interleaved', 'paired', or 'single'
                'seq_tech': sequencing_tech,
                'reads_ref': KBase object ref for downstream convenience,
                'reads_name': KBase object name for downstream convenience,
                'rev_file': path_to_fastq_file, #only if paired end
        }
        """
        rds_params = copy.deepcopy(input_params)
        wsname = rds_params[self.PARAM_IN_WS]
        rds_refs = []
        rds_data = []

        # jump_libraries grouped params
        if 'jump_libraries' in rds_params and rds_params['jump_libraries']:
            self._unique_prefix_check('jp_prefix', rds_params['jump_libraries'])
            for rds_lib in rds_params['jump_libraries']:
                if 'jp_id' in rds_lib:
                    rds_refs.append(rds_lib['jp_id'])
            rds_data = self._getKBReadsInfo(wsname, rds_refs)

            for rds_lib in rds_params['jump_libraries']:
                i = 0
                for rds in rds_data:
                    i += 1
                    if ('jp_id' in rds_lib and rds_lib['jp_id'] == rds['reads_ref']):
                        if 'jp_prefix' in rds_lib:
                            rds['jp_prefix'] = rds_lib['jp_prefix'][:1] + str(i)
                        else:
                            raise ValueError("Parameter jp_prefix is required for reads {}".format(
                                rds[reads_ref]))
                        if 'jp_mean' in rds_lib:
                            rds['jp_mean'] = rds_lib['jp_mean']
                        else:
                            raise ValueError(
                                "Parameter jp_mean is required for reads {}".format(rds[reads_ref]))
                        if 'jp_stdev' in rds_lib:
                            rds['jp_stdev'] = rds_lib['jp_stdev']
                        else:
                            raise ValueError("Parameter pe_stdev is required for reads {}".format(
                                rds[reads_ref]))
        return rds_data

    def _getKBReadsInfo(self, wsname, reads_refs):
        """
        _getKBReadsInfo--from a set of given KBase reads refs, fetches the corresponding reads info
        with as deinterleaved fastq files and returns a list of reads data in the following structure:
        reads_data = {
                'fwd_file': path_to_fastq_file,
                'type': reads_type, #('interleaved', 'paired', or 'single'
                'seq_tech': sequencing_tech,
                'reads_ref': KBase object ref for downstream convenience,
                'reads_name': KBase object name for downstream convenience,
                'rev_file': path_to_fastq_file, #only if paired end
        }
        """
        obj_ids = []
        for r in reads_refs:
            if r:
                obj_ids.append({'ref': r if '/' in r else (wsname + '/' + r)})

        if not obj_ids:
            return []

        ws_info = self.ws_client.get_object_info_new({'objects': obj_ids})
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
            reads = self.ru.download_reads({
                        'read_libraries': reads_params,
                        'interleaved': 'false'
                        })['files']
        except ServerError as se:
            log('logging stacktrace from dynamic client error')
            log(se.data)
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

        # log('Downloaded reads data from KBase:\n' + pformat(reads))
        reads_data = []
        for ref in reads_refs:
            reads_name = reftoname[ref]
            f = reads[ref]['files']
            seq_tech = reads[ref]['sequencing_tech']
            rds_info = {
                'fwd_file': f['fwd'],
                'reads_ref': ref,
                'type': f['type'],
                'seq_tech': seq_tech,
                'reads_name': reads_name
            }
            if f.get('rev', None) is not None:
                rds_info['rev_file'] = f['rev']
            reads_data.append(rds_info)

        return reads_data

    def _generate_output_file_list(self, out_dir):
        """
        _generate_output_file_list: zip result files and generate file_links for report
        """
        log('start packing result files')

        output_files = list()

        output_directory = os.path.join(self.proj_dir, str(uuid.uuid4()))
        _mkdir_p(output_directory)
        masurca_output = os.path.join(output_directory, 'masurca_output.zip')
        self._zip_folder(out_dir, masurca_output)

        output_files.append({'path': masurca_output,
                             'name': os.path.basename(masurca_output),
                             'label': os.path.basename(masurca_output),
                             'description': 'Output file(s) generated by MaSuRCA'})

        return output_files

    def _zip_folder(self, folder_path, output_path):
        """
        _zip_folder: Zip the contents of an entire folder (with that folder included
        in the archive). Empty subfolders could be included in the archive as well
        if the commented portion is used.
        """
        with zipfile.ZipFile(output_path, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as ziph:
            for root, folders, files in os.walk(folder_path):
                for f in files:
                    absolute_path = os.path.join(root, f)
                    relative_path = os.path.join(os.path.basename(root), f)
                    # print "Adding {} to archive.".format(absolute_path)
                    ziph.write(absolute_path, relative_path)

        print("{} created successfully.".format(output_path))
        # with zipfile.ZipFile(output_path, "r") as f:
        #    print 'Checking the zipped file......\n'
        #    for info in f.infolist():
        #        print info.filename, info.date_time, info.file_size, info.compress_size

    def _load_stats(self, input_file_name):
        log('Starting conversion of FASTA to KBaseGenomeAnnotations.Assembly')
        log('Building Object.')
        if not os.path.isfile(input_file_name):
            raise Exception('The input file name {0} is not a file!'.format(input_file_name))
        with open(input_file_name, 'r') as input_file_handle:
            contig_id = None
            sequence_len = 0
            fasta_dict = dict()
            first_header_found = False
            # Pattern for replacing white space
            pattern = re.compile(r'\s+')
            for current_line in input_file_handle:
                if (current_line[0] == '>'):
                    # found a header line
                    # Wrap up previous fasta sequence
                    if not first_header_found:
                        first_header_found = True
                    else:
                        fasta_dict[contig_id] = sequence_len
                        sequence_len = 0
                    fasta_header = current_line.replace('>', '').strip()
                    try:
                        contig_id = fasta_header.strip().split(' ', 1)[0]
                    except:
                        contig_id = fasta_header.strip()
                else:
                    sequence_len += len(re.sub(pattern, '', current_line))
        # wrap up last fasta sequence
        if not first_header_found:
            raise Exception("There are no contigs in this file")
        else:
            fasta_dict[contig_id] = sequence_len
        return fasta_dict

    def _valid_string(self, s_str, is_ref=False):
        is_valid = isinstance(s_str, str) and len(s_str.strip()) > 0
        if is_valid and is_ref:
            is_valid = self._check_reference(s_str)
        return is_valid

    def _check_reference(self, ref):
        """
        Tests the given ref string to make sure it conforms to the expected
        object reference format. Returns True if it passes, False otherwise.
        """
        obj_ref_regex = re.compile("^(?P<wsid>\d+)\/(?P<objid>\d+)(\/(?P<ver>\d+))?$")
        ref_path = ref.strip().split(";")
        for step in ref_path:
            if not obj_ref_regex.match(step):
                return False
        return True

    def _check_ref_type(self, ref, allowed_types):
        """
        Validates the object type of ref against the list of allowed types. If it passes, this
        returns True, otherwise False.
        Really, all this does is verify that at least one of the strings in allowed_types is
        a substring of the ref object type name.
        Ex1:
        ref = "KBaseGenomes.Genome-4.0"
        allowed_types = ["assembly", "KBaseFile.Assembly"]
        returns False
        Ex2:
        ref = "KBaseGenomes.Genome-4.0"
        allowed_types = ["assembly", "genome"]
        returns True
        """
        obj_type = self._get_object_type(ref).lower()
        for t in allowed_types:
            if t.lower() in obj_type:
                return True
        return False

    def _get_object_type(self, ref):
        """
        Fetches and returns the typed object name of ref from the given workspace url.
        If that object doesn't exist, or there's another Workspace error, this raises a
        RuntimeError exception.
        """
        info = self.ws_client.get_object_info3({'objects': [{'ref': ref}]})
        obj_info = info.get('infos', [[]])[0]
        if len(obj_info) == 0:
            raise RuntimeError("An error occurred while fetching type info from the Workspace. "
                               "No information returned for reference {}".format(ref))
        return obj_info[2]

    def _get_fasta_from_assembly(self, assembly_ref):
        """
        From an assembly or contigset, this uses a data file to build a FASTA file
        and return the path to it.
        """
        allowed_types = ['KBaseFile.Assembly',
                         'KBaseGenomeAnnotations.Assembly',
                         'KBaseGenomes.ContigSet']
        if not self._check_ref_type(assembly_ref, allowed_types):
            raise ValueError(
                "The reference {} cannot be used to fetch a FASTA file".format(assembly_ref))
        au = AssemblyUtil(self.callback_url)
        return au.get_assembly_as_fasta({'ref': assembly_ref})

    def generate_report(self, contig_file_name, params, out_dir, wsname):
        log('Generating and saving report')

        contig_file_with_path = os.path.join(out_dir, contig_file_name)
        fasta_stats = self._load_stats(contig_file_with_path)
        lengths = [fasta_stats[contig_id] for contig_id in fasta_stats]

        assembly_ref = params[self.PARAM_IN_WS] + '/' + params[self.PARAM_IN_CS_NAME]

        report_text = ''
        report_text += 'MaSuRCA results saved to: ' + wsname + '/' + out_dir + '\n'
        report_text += 'Assembly saved to: ' + assembly_ref + '\n'
        report_text += 'Assembled into ' + str(len(lengths)) + ' contigs.\n'
        report_text += 'Avg Length: ' + str(sum(lengths) / float(len(lengths))) + ' bp.\n'

        # compute a simple contig length distribution
        bins = 10
        counts, edges = np.histogram(lengths, bins)
        report_text += 'Contig Length Distribution (# of contigs -- min to max ' + 'basepairs):\n'
        for c in range(bins):
            report_text += ('   ' + str(counts[c]) + '\t--\t' + str(edges[c]) + ' to ' +
                            str(edges[c + 1]) + ' bp\n')
        print('Running QUAST')
        quastret = self.kbq.run_QUAST(
            {'files': [{'path': contig_file_with_path, 'label': params[self.PARAM_IN_CS_NAME]}]})

        output_files = self._generate_output_file_list(out_dir)

        print('Saving report')
        report_output = self.kbr.create_extended_report(
            {'message': report_text,
             'objects_created': [{'ref': assembly_ref, 'description': 'Assembled contigs'}],
             'direct_html_link_index': 0,
             'file_links': output_files,
             'html_links': [{'shock_id': quastret['shock_id'],
                             'name': 'report.html',
                             'label': 'QUAST report'}
                            ],
             'report_object_name': 'kb_masurca_report_' + str(uuid.uuid4()),
             'workspace_name': params[self.PARAM_IN_WS]})
        report_name = report_output['name']
        report_ref = report_output['ref']
        return report_name, report_ref

    def validate_params(self, params):
        """
        validate_params: checks params passed to run_masurca_app method and set default values
        """
        # log('Start validating run_masurca_app parameters:\n{}'.format(json.dumps(params, indent=1)))

        # check for mandatory parameters
        if params.get(self.PARAM_IN_WS, None) is None:
            raise ValueError(self.PARAM_IN_WS + ' parameter is mandatory')
        if self.PARAM_IN_THREADN not in params:
            raise ValueError(self.PARAM_IN_THREADN + ' parameter is mandatory')

        if params.get(self.PARAM_IN_JF_SIZE, None) is None:
            raise ValueError(self.PARAM_IN_JF_SIZE + ' parameter is mandatory')
        if params.get(self.PARAM_IN_READS_LIBS, None) is None:
            raise ValueError(self.PARAM_IN_READS_LIBS + ' parameter is mandatory')
        if type(params[self.PARAM_IN_READS_LIBS]) != list:
            raise ValueError(self.PARAM_IN_READS_LIBS + ' must be a list')

        if (params.get(self.PARAM_IN_CS_NAME, None) is None or
                not self._valid_string(params[self.PARAM_IN_CS_NAME])):
            raise ValueError("Parameter output_contigset_name is required and " +
                             "must be a valid Workspace object string, " +
                             "not {}".format(params.get(self.PARAM_IN_CS_NAME, None)))
        if ('pe_prefix' not in params):
            params['pe_prefix'] = 'pe'
        if ('pe_mean' not in params or type(params['pe_mean']) != int):
            params['pe_mean'] = 180
        if ('pe_stdev' not in params or type(params['pe_stdev']) != int):
            params['pe_stdev'] = 20
        if ('jp_prefix' not in params):
            params['jp_prefix'] = 'sh'
        if ('dna_source' in params):
            dna_src = params.get('dna_source')
            if dna_src == 'bacteria':
                params['limit_jump_coverage'] = 60
                params['cgwErrorRate'] = 0.25
            else:
                params['limit_jump_coverage'] = 300
                params['cgwErrorRate'] = 0.15

        if params.get('create_report', None) is None:
            params['create_report'] = 0

        return params

    def construct_masurca_assembler_cfg(self, params):
        # STEP 1: get the working folder housing the config.txt file and the masurca results
        wsname = params[self.PARAM_IN_WS]
        config_file_path = os.path.join(self.proj_dir, 'config.txt')

        # STEP 2.1: retrieve the reads data from input parameter
        pe_reads_data = self._getReadsInfo_PE(params)
        jp_reads_data = []
        if params.get(self.PARAM_IN_JUMP_LIBS, None):
            jp_reads_data = self._getReadsInfo_JP(params)
            if ('jp_mean' not in params or type(params['jp_mean']) != int):
                params['jp_mean'] = 3600
            if ('jp_stdev' not in params or type(params['jp_stdev']) != int):
                params['jp_stdev'] = 200

        # STEP 2.2: PACBIO reads must be in a single FASTA file and supplied as PACBIO=reads.fa;
        assbl_types = ['KBaseFile.Assembly',
                       'KBaseGenomeAnnotations.Assembly',
                       'KBaseGenomes.ContigSet']
        reads_types = ['KBaseAssembly.SingleEndLibrary',
                       'KBaseFile.SingleEndLibrary',
                       'KBaseAssembly.PairedEndLibrary',
                       'KBaseFile.PairedEndLibrary']
        pb_reads_file = ''
        if params.get('pacbio_reads', None):
            pb_ref = params['pacbio_reads']
            if self._check_ref_type(pb_ref, assbl_types):
                pb_reads_file = (self._get_fasta_from_assembly(pb_ref)).get('path', '')
            else:
                if self._check_ref_type(pb_ref, reads_types):
                    pb_rd = self._getKBReadsInfo(wsname, [pb_ref])
                    pb_reads_file = pb_rd[0]['fwd_file']
                    if pb_rd[0].get('rev_file', None):
                        pb_reads_file += ' ' + pb_rd[0]['rev_file']

        # STEP 2.3: NANOPORE reads must be in a single FASTA/FASTQ file and supplied as NANOPORE=reads.fa
        np_reads_file = ''
        if params.get('nanopore_reads', None):
            np_ref = params['nanopore_reads']
            if self._check_ref_type(np_ref, assbl_types):
                np_reads_file = (self._get_fasta_from_assembly(np_ref)).get('path', '')
            else:
                if self._check_ref_type(np_ref, reads_types):
                    np_rd = self._getKBReadsInfo(wsname, [np_ref])
                    np_reads_file = np_rd[0]['fwd_file']
                    if np_rd[0].get('rev_file', None):
                        np_reads_file += ' ' + np_rd[0]['rev_file']

        # STEP 2.4: any OTHER sequence data (454, Sanger, Ion torrent, etc) must be first
        # converted into Celera Assembler compatible .frg files
        # (see http://wgsassembler.sourceforge.com) and supplied as OTHER=file.frg
        other_frg = ''
        if params.get('other_frg_file', None):
            other_frg = params['other_frg_file']

        # STEP 3: construct and save the config.txt file for running masurca
        try:
            # STEP 3.1: replace the 'DATA...END' portion of the config_template.txt file
            data_str = self._get_data_portion(pe_reads_data, jp_reads_data,
                                              pb_reads_file, np_reads_file, other_frg)
            if data_str == '':  # no reads libraries are specified, no further actions
                return ''

            config_template = ''
            with codecs.open(
                             os.path.join(os.path.dirname(__file__), 'config_template.txt'),
                             mode='r', encoding='utf-8') as config_template_file:
                config_template = config_template_file.read()

            begin_patn1 = "DATA\n"
            end_patn1 = "END\nPARAMETERS\n"
            config_with_data = self._replaceSectionText(config_template, begin_patn1,
                                                        end_patn1, data_str)
            # log("\n***After DATA section replacement:\n{}\nSaved at {}".format(
            #             config_with_data.encode('utf-8').decode('utf-8'), config_file_path))

            with codecs.open(config_file_path, mode='w', encoding='utf-8') as config_file:
                config_file.write(config_with_data)

            # STEP 3.2: replace the 'PARAMETERS...END' portion of the config_file file saved in STEP 3.1
            param_str = self._get_parameters_portion(params)
            if param_str == '':  # no parameters are specified, no further actions
                return ''

            previous_config = ''
            with codecs.open(config_file_path, mode='r', encoding='utf-8') as previous_config_file:
                previous_config = previous_config_file.read()

            begin_patn2 = "PARAMETERS\n"
            end_patn2 = "END\n"
            final_config = self._replaceSectionText(previous_config, begin_patn2,
                                                    end_patn2, param_str)
            log("\n***Configuration file content:\n{}\nSaved at {}".format(
                        final_config.encode('utf-8').decode('utf-8'), config_file_path))

            with codecs.open(config_file_path, mode='w', encoding='utf-8') as config_file:
                config_file.write(final_config)
        except IOError as ioerr:
            log('Creation of the config.txt file raised error:\n')
            pprint(ioerr)
            return ''
        else:
            return config_file_path

    def generate_assemble_script(self, config_file):
        exit_code = 1
        if os.path.isfile(config_file):
            f_dir, f_nm = os.path.split(config_file)
            m_cmd = [self.MaSuRCA_BIN]
            m_cmd.append(config_file)
            exit_code = self.prog_runner.run(m_cmd, f_dir)

            if exit_code == 0:
                log('Created the assemble.sh file at {}.\n'.format(
                    os.path.join(f_dir, 'assemble.sh')))
                return os.path.join(f_dir, 'assemble.sh')
            else:
                log('Creation of the assemble.sh file failed.\n')
                return ''
        else:
            log("The config file {} is not found.\n".format(config_file))
            log('NO assemble.sh file created.\n')
        return ''

    def run_assemble(self, asmbl_file):
        exit_code = 1
        if os.path.isfile(asmbl_file):
            log("The assemble.sh file exists at {}\n".format(asmbl_file))
            f_dir, f_nm = os.path.split(asmbl_file)
            a_cmd = ['/bin/bash']
            a_cmd.append(asmbl_file)
            log("The working directory is {}\n".format(f_dir))
            log("The assembling command is {}\n".format(' '.join(a_cmd)))

            p = subprocess.Popen(a_cmd, cwd=f_dir, shell=False)
            exit_code = p.wait()
            log('Return code: ' + str(exit_code))

            if p.returncode != 0:
                raise ValueError('Error running assemble.sh, return code: ' +
                                 str(p.returncode) + '\n')
            else:
                exit_code = p.returncode
        else:
            log("The assemble.sh file {} is not found.".format(asmbl_file))
        return exit_code

    def save_assembly(self, contig_fa, wsname, a_name):
        if os.path.isfile(contig_fa):
            log('Uploading FASTA file to Assembly...')
            self.au.save_assembly_from_fasta(
                            {'file': {'path': contig_fa},
                             'workspace_name': wsname,
                             'assembly_name': a_name})
        else:
            log("The contig file {} is not found.".format(contig_fa))
