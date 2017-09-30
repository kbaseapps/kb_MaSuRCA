# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import requests
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
from ReadsUtils.baseclient import ServerError
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil

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
        assembly_ref = assemblyUtil.save_assembly_from_fasta({'file': {'path': filename},
                                                              'workspace_name': self.getWsName(),
                                                              'assembly_name': obj_name
                                                              })
        return assembly_ref

    # borrowed from Megahit - call this method to get the WS object info of a Paired End Library (will
    # upload the example data if this is the first time the method is called during tests)
    def loadPairedEndReads(self):
        if hasattr(self.__class__, 'pairedEndLibInfo'):
            return self.__class__.pairedEndLibInfo
        # 1) upload files to shock
        shared_dir = "/kb/module/work/tmp"
        forward_data_file = '../test/testReads/small.forward.fq'
        forward_file = os.path.join(shared_dir, os.path.basename(forward_data_file))
        shutil.copy(forward_data_file, forward_file)
        reverse_data_file = '../test/testReads/small.reverse.fq'
        reverse_file = os.path.join(shared_dir, os.path.basename(reverse_data_file))
        shutil.copy(reverse_data_file, reverse_file)

        ru = ReadsUtils(os.environ['SDK_CALLBACK_URL'])
        pe_reads_ref = ru.upload_reads({'fwd_file': forward_file, 'rev_file': reverse_file,
                                          'sequencing_tech': 'artificial reads',
                                          'interleaved': 0, 'wsname': self.getWsName(),
                                          'name': 'test_pe_reads'})['obj_ref']

        self.__class__.pe_reads_ref = pe_reads_ref
        print('Loaded PairedEndReads: ' + pe_reads_ref)
        new_obj_info = self.wsClient.get_object_info_new({'objects': [{'ref': pe_reads_ref}]})
        self.__class__.pairedEndLibInfo = new_obj_info[0]
        pprint (pformat(new_obj_info))
        #return new_obj_info[0]
        return pe_reads_ref

    def loadSEReads(self, reads_file_path):
        #if hasattr(self.__class__, 'reads_ref'):
            #return self.__class__.reads_ref
        se_reads_name = os.path.basename(reads_file_path)
        fq_path = os.path.join(self.scratch, se_reads_name)
        shutil.copy(reads_file_path, fq_path)

        ru = ReadsUtils(self.callback_url)
        reads_ref = ru.upload_reads({'fwd_file': fq_path,
                                        'wsname': self.getWsName(),
                                        'name': se_reads_name.split('.')[0],
                                        'sequencing_tech': 'rnaseq reads'})['obj_ref']
        #self.__class__.reads_ref = reads_ref
        return reads_ref



    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    # Uncomment to skip this test
    # @unittest.skip("skipped test_run_masurca")
    def test_run_masurca(self):
        # First load a test FASTA file as an KBase Assembly
        se_lib_ref = self.loadSEReads(os.path.join('../test/testReads', 'small.forward.fq'))
        pe_lib_ref = self.loadPairedEndReads()
        m_params = {
                'workspace_name': self.getWsName(),
                'reads_libraries': [pe_lib_ref],
                "jf_size": 100000000,
                "pe_prefix": "pe",
                "pe_mean": 180,
                "pe_stdev": 20,
                #"jp_prefix": "sh",
                #"jp_mean": 3600,
                #"jp_stdev": 200,
                'jump_libraries': None,
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
                'create_report': 1
        }
        # Second, call your implementation
        ret = self.getImpl().run_masurca(self.getContext(), m_params)

        # Validate the returned data
        #self.assertEqual(ret[0]['n_initial_contigs'], 3)
        #self.assertEqual(ret[0]['n_contigs_removed'], 1)
        #self.assertEqual(ret[0]['n_contigs_remaining'], 2)

    @unittest.skip("skipped test_run_masurca_err1")
    def test_run_masurca_err1(self):
        with self.assertRaises(ValueError) as errorContext:
            self.getImpl().run_masurca(self.getContext(),
                                          {'workspace_name': self.getWsName(),
                                           'assembly_input_ref': '1/fake/3',
                                           'min_length': '-10'})
        self.assertIn('min_length parameter cannot be negative', str(errorContext.exception))

    @unittest.skip("skipped test_run_masurca_err2")
    def test_run_masurca_err2(self):
        with self.assertRaises(ValueError) as errorContext:
            self.getImpl().run_masurca(self.getContext(),
                                          {'workspace_name': self.getWsName(),
                                           'assembly_input_ref': '1/fake/3',
                                           'min_length': 'ten'})
        self.assertIn('Cannot parse integer from min_length parameter', str(errorContext.exception))
