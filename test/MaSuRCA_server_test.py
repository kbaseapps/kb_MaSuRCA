# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import shutil

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint, pformat  # noqa: F401

from biokbase.workspace.client import Workspace as workspaceService
from MaSuRCA.MaSuRCAImpl import kb_MaSuRCA
from MaSuRCA.MaSuRCAServer import MethodContext
from MaSuRCA.authclient import KBaseAuth as _KBaseAuth
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.AssemblyUtilClient import AssemblyUtil

from MaSuRCA.core.masurca_assembler import MaSuRCA_Assembler
from MaSuRCA.core.masurca_utils import masurca_utils


class MaSuRCATest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_MaSuRCA'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'MaSuRCA',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = kb_MaSuRCA(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

        cls.masurca_PROJECT_DIR = 'masurca_outputs'
        cls.scratch = cls.cfg['scratch']
        if not os.path.exists(cls.scratch):
            os.makedirs(cls.scratch)
        cls.masurca_prjdir = os.path.join(cls.scratch, cls.masurca_PROJECT_DIR)
        cls.masurca_assembler = MaSuRCA_Assembler(cls.cfg, cls.ctx.provenance)
        cls.masurca_utils = masurca_utils(cls.masurca_prjdir, cls.cfg)

    @classmethod
    def make_ref(self, object_info):
        return str(object_info[6]) + '/' + str(object_info[0]) + \
            '/' + str(object_info[4])

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_MaSuRCA_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def load_fasta_file(self, filename, obj_name, contents):
        f = open(filename, 'w')
        f.write(contents)
        f.close()
        assemblyUtil = AssemblyUtil(self.callback_url)
        assembly_ref = assemblyUtil.save_assembly_from_fasta(
            {'file': {'path': filename},
             'workspace_name': self.getWsName(),
             'assembly_name': obj_name})
        return assembly_ref

    def loadAssembly(self, fa_file_path):
        # if hasattr(self.__class__, 'assembly_ref'):
        #   return self.__class__.assembly_ref

        assembly_nm = os.path.basename(fa_file_path)
        fasta_path = os.path.join(self.scratch, assembly_nm)
        shutil.copy(fa_file_path, fasta_path)
        au = AssemblyUtil(self.callback_url)
        assembly_ref = au.save_assembly_from_fasta(
            {'file': {'path': fasta_path},
             'workspace_name': self.getWsName(),
             'assembly_name': assembly_nm.split('.')[0]})
        # self.__class__.assembly_ref = assembly_ref
        print('Loaded Assembly:{} with ref of{}.'.format(
            assembly_nm, assembly_ref))
        return assembly_ref

    def loadSEReads(self, reads_file_path):
        #  if hasattr(self.__class__, 'reads_ref'):
        #      return self.__class__.reads_ref
        se_reads_name = os.path.basename(reads_file_path)
        fq_path = os.path.join(self.scratch, se_reads_name)
        shutil.copy(reads_file_path, fq_path)

        ru = ReadsUtils(self.callback_url)
        reads_ref = ru.upload_reads({'fwd_file': fq_path,
                                     'wsname': self.getWsName(),
                                     'name': se_reads_name.split('.')[0],
                                     'sequencing_tech': 'kb reads'})['obj_ref']
        #  self.__class__.reads_ref = reads_ref
        return reads_ref

    #  borrowed from Megahit - call this method to get the WS object info of a Paired End Library
    # (will upload the example data if this is the first time the method is called during tests)
    def loadPairedEndReads(self, forward_data_file, reverse_data_file, reads_name):
        # if hasattr(self.__class__, 'pairedEndLibInfo'):
        #    return self.__class__.pairedEndLibInfo
        # 1) upload files to shock
        # shared_dir = "/kb/module/work/tmp"
        # forward_data_file = '../test/testReads/small.forward.fq'
        forward_file = os.path.join(
            self.scratch, os.path.basename(forward_data_file))
        shutil.copy(forward_data_file, forward_file)
        # reverse_data_file = '../test/testReads/small.reverse.fq'
        reverse_file = os.path.join(
            self.scratch, os.path.basename(reverse_data_file))
        shutil.copy(reverse_data_file, reverse_file)

        ru = ReadsUtils(os.environ['SDK_CALLBACK_URL'])
        pe_reads_ref = ru.upload_reads(
            {'fwd_file': forward_file, 'rev_file': reverse_file,
             'sequencing_tech': 'artificial reads',
             'interleaved': 0, 'wsname': self.getWsName(),
             'name': reads_name})['obj_ref']

        # self.__class__.pe_reads_ref = pe_reads_ref
        print('Loaded PairedEndReads: ' + pe_reads_ref)
        new_obj_info = self.wsClient.get_object_info_new(
            {'objects': [{'ref': pe_reads_ref}]})
        self.__class__.pairedEndLibInfo = new_obj_info[0]
        pprint(pformat(new_obj_info))
        # return new_obj_info[0]
        return pe_reads_ref

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    # Uncomment to skip this test
    # @unittest.skip("skipped test_run_masurca_assemblerr_with_longreads")
    def test_run_masurca_assembler_with_longreads(self):
        # First load a test FASTA file as an KBase Assembly
        # se_lib_ref1 = self.loadSEReads(
        # os.path.join('../test/testReads', 'short_reads_1.fastq'))
        # se_lib_ref2 = self.loadSEReads(
        # os.path.join('../test/testReads', 'short_reads_2.fastq'))
        se_lib_ref3 = self.loadSEReads(
            os.path.join('../test/testReads', 'testreads.wt.fastq'))
        asmbl_ref1 = self.loadAssembly(
            os.path.join('../test/testReads', 'test_reference.fa'))
        asmbl_ref2 = self.loadAssembly(
            os.path.join('../test/testReads', 'reference.fasta'))
        pe_lib_ref1 = self.loadPairedEndReads(
            '../test/testReads/small.forward.fq',
            '../test/testReads/small.reverse.fq', 'small')
        pe_lib_ref2 = self.loadPairedEndReads(
            '../test/testReads/short_reads_1.fastq',
            '../test/testReads/short_reads_2.fastq', 'short')

        m_params = {
            "workspace_name": self.getWsName(),
            "jf_size": 100000000,
            "reads_libraries": [{
                "pe_id": pe_lib_ref1,
                "pe_prefix": "pe",
                "pe_mean": 180,
                "pe_stdev": 20
            }, {
                "pe_id": se_lib_ref3,
                "pe_prefix": "se",
                "pe_mean": 180,
                "pe_stdev": 20
            }],
            "jump_libraries": [{
                "jp_id": pe_lib_ref2,
                "jp_prefix": "sh",
                "jp_mean": 3600,
                "jp_stdev": 200
            }],
            "pacbio_reads": asmbl_ref1,
            "nanopore_reads": asmbl_ref2,
            "other_frg_file": "",
            "output_contigset_name": "withLongReads_masurca.contigs",
            "graph_kmer_size": None,
            "use_linking_mates": 1,
            "limit_jump_coverage": 60,
            "cgwErrorRate": 0.25,
            "close_gaps": 0,
            "soap_assembly": 0,
            "do_homopolymer_trim": 0,
            "kmer_count_threshold": 1,
            "num_threads": 2,
            "create_report": 1
        }
        # Because the selected reads files are random,
        # masurca reported it 'Could not find file final.genome.scf.fasta!'
        with self.assertRaises(ValueError) as errorContext:
            self.getImpl().run_masurca_assembler(self.getContext(), m_params)
            self.assertIn('masurca assemble process failed',
                          str(errorContext.exception))

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    # Uncomment to skip this test
    # @unittest.skip("skipped test_run_masurca_assembler_no_longreads")
    def test_run_masurca_assembler_no_longreads(self):
        # First load a test FASTA file as an KBase Assembly
        se_lib_ref = self.loadSEReads(
            os.path.join('../test/testReads', 'testreads.wt.fastq'))
        pe_lib_ref1 = self.loadPairedEndReads(
            '../test/testReads/small.forward.fq',
            '../test/testReads/small.reverse.fq', 'small')
        pe_lib_ref2 = self.loadPairedEndReads(
            '../test/testReads/short_reads_1.fastq',
            '../test/testReads/short_reads_2.fastq', 'short')
        m_params = {
            'workspace_name': self.getWsName(),
            'reads_libraries': [{
                "pe_id": pe_lib_ref1,
                "pe_prefix": "p1",
                "pe_mean": 180,
                "pe_stdev": 20
            }, {
                "pe_id": pe_lib_ref2,
                "pe_prefix": "p2",
                "pe_mean": 180,
                "pe_stdev": 20
            }],
            'jump_libraries': [{
                "jp_id": se_lib_ref,
                "jp_prefix": "sh",
                "jp_mean": 3600,
                "jp_stdev": 200
            }],
            "jf_size": 100000000,
            "pacbio_reads": None,
            "other_frg_file": None,
            "output_contigset_name": "masurca.contigs",
            "graph_kmer_size": "auto",
            "use_linking_mates": 1,
            "limit_jump_coverage": 60,
            "cgwErrorRate": 0.25,
            "close_gaps": 0,
            "soap_assembly": 0,
            "do_homopolymer_trim": 0,
            "kmer_count_threshold": 1,
            'num_threads': 2,
            'create_report': 0
        }
        # Because the selected reads files are random,
        # masurca reported it 'Could not find file final.genome.scf.fasta!'
        with self.assertRaises(ValueError) as errorContext:
            self.getImpl().run_masurca_assembler(self.getContext(), m_params)
            self.assertIn('Failed to generate assemble.sh file!',
                          str(errorContext.exception))

    # @unittest.skip("skipped test_masurca_utils_validate_params")
    def test_masurca_utils_validate_params(self):
        with self.assertRaises(ValueError) as errorContext:
            self.masurca_utils.validate_params(
                {'workspace_name': self.getWsName(),
                 'jf_size': 100000000,
                 'num_threads': 2,
                 'output_contigset_name': 'masurca.contigs'})
        self.assertIn('reads_libraries parameter is mandatory',
                      str(errorContext.exception))
        with self.assertRaises(ValueError) as errorContext:
            self.masurca_utils.validate_params(
                {'workspace_name': self.getWsName(),
                 'reads_libraries': [{'pe_id': '1/fake/3'}],
                 'jf_size': 100000000,
                 'output_contigset_name': 'masurca.contigs'})
        self.assertIn('num_threads parameter is mandatory',
                      str(errorContext.exception))
        with self.assertRaises(ValueError) as errorContext:
            self.masurca_utils.validate_params(
                {'workspace_name': self.getWsName(),
                 'reads_libraries': [{'pe_id': '1/fake/3'}],
                 'num_threads': 2,
                 'output_contigset_name': 'masurca.contigs'})
        self.assertIn('jf_size parameter is mandatory',
                      str(errorContext.exception))
        with self.assertRaises(ValueError) as errorContext:
            self.masurca_utils.validate_params(
                {'workspace_name': self.getWsName(),
                 'reads_libraries': [{'pe_id': '1/fake/3'}],
                 'jf_size': 100000000,
                 'num_threads': 2})
        self.assertIn('output_contigset_name is mandatory',
                      str(errorContext.exception))
        with self.assertRaises(ValueError) as errorContext:
            self.masurca_utils.validate_params(
                {'workspace_name': self.getWsName(),
                 'reads_libraries': {'pe_id': '1/fake/3'},
                 'jf_size': 100000000,
                 'num_threads': 2,
                 'output_contigset_name': 'masurca.contigs',
                 'min_length': 100})
            self.assertIn('must be a list',
                          str(errorContext.exception))
