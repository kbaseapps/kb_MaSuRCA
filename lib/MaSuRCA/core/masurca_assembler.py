import os
import re
import time
import uuid
import json
import shutil
import subprocess
import numpy as np
from datetime import datetime
from Bio import SeqIO
from pprint import pprint, pformat
import traceback
import zipfile
import multiprocessing
import psutil

#from DataFileUtil.DataFileUtilClient import DataFileUtil
#from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil

from KBaseReport.KBaseReportClient import KBaseReport
from KBaseReport.baseclient import ServerError as _RepError
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from kb_quast.kb_quastClient import kb_quast
from kb_quast.baseclient import ServerError as QUASTError
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from ReadsUtils.baseclient import ServerError
from Workspace.WorkspaceClient import Workspace as workspace

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class MaSuRCA_Assembler(object):
    INVALID_WS_OBJ_NAME_RE = re.compile('[^\\w\\|._-]')
    INVALID_WS_NAME_RE = re.compile('[^\\w:._-]')

    PARAM_IN_CS_NAME = 'output_contigset_name'
    PARAM_IN_READ_LIB = 'read_libraries'
    PARAM_IN_JUMP_LIB = 'jump_libraries'
    PARAM_IN_PACBIOS = 'pacbios'
    PARAM_IN_OTHER = 'read_libraries'
    MaSuRCAR_PROJECT_DIR = 'masurca_project_dir'

    def __init__(self, config, provenance):
        """
        __init__
        """
        # BEGIN_CONSTRUCTOR
        self.ws_url = config["workspace-url"]
        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.shock_url = config['shock-url']
        self.dfu = DataFileUtil(self.callback_url, service_ver='beta')
        self.gfu = GenomeFileUtil(self.callback_url)
        self.au = AssemblyUtil(self.callback_url)
        self.ws = Workspace(self.ws_url, token=self.token)
        self.provenance = provenance

        self.scratch = os.path.join(config['scratch'], str(uuid.uuid4()))
        self._mkdir_p(self.scratch)

        self.masurca_version = 'MaSuRCA-' + os.environ['M_VERSION']

        # from the provenance, extract out the version to run by exact hash if possible
        self.my_version = 'release'
        if len(provenance) > 0:
            if 'subactions' in provenance[0]:
                self.my_version = self.get_version_from_subactions('kb_MaSuRCA', provenance[0]['subactions'])
        print('Running kb_MaSuRCA version = ' + self.my_version)
        # END_CONSTRUCTOR
        pass


    def run_masurca(self, params):
        # 0. create the masurca project folder
        if self.proj_dir is None:
            m_dir = self.create_star_dirs(self.scratch)
            self.proj_dir = m_dir
        wsname = params['workspace_name']

        cpus = min(params.get('num_threads'), psutil.cpu_count())

        # 1. validate & process the input parameters
        validated_params = self.process_params(params)

        # 2. convert the input parameters (from refs to file paths, especially)
        input_params = self.convert_params(validated_params)

        # 3. create the configuratio file 
        config_file = self.create_config_file(input_params)

        # 4. run masurca against the configuratio file to generate the assemble.sh script 
        assemble_file = self.generate_assemble_script(config_file)

        # 5. run the assemble.sh script to do the heavy-lifting
        self.run_assemble(assemble_file)

        # 6. save the assembly to KBase if everything has gone well
        self.log('Uploading FASTA file to Assembly')
        assemblyUtil = AssemblyUtil(self.callbackURL, token=ctx['token'], service_ver='release')
        output_contigs = os.path.join(self.proj_dir, 'name_of_fa_file')
        assemblyUtil.save_assembly_from_fasta(
                        {'file': {'path': output_contigs},
                        'workspace_name': wsname,
                        'assembly_name': params[self.PARAM_IN_CS_NAME]
                        })

        # 7. report the final results
        returnVal = {
            "report_ref": None,
            "report_name": None
        }

        report_name, report_ref = self.generate_report(output_contigs, params, proj_dir, wsname)
        returnVal = {'report_name': report_name, 'report_ref': report_ref}

        return returnVal


    def create_config_file(self, params):
        '''creating the configuration file for MaSuRCA to generate the assemble.sh script'''
        config_file_with_path = ''

        return config_file_with_path


    def create_proj_dir(self, home_dir):
        '''creating the project directory for MaSuRCA'''
        prjdir = os.path.join(home_dir, self.MaSuRCAR_PROJECT_DIR)
        self._mkdir_p(prjdir)

        return prjdir


    def get_version_from_subactions(self, module_name, subactions):
        # go through each sub action looking for
        if not subactions:
            return 'dev' #'release'  # default to release if we can't find anything
        for sa in subactions:
            if 'name' in sa:
                if sa['name'] == module_name:
                    # local-docker-image implies that we are running in kb-test, so return 'dev'
                    if sa['commit'] == 'local-docker-image':
                        return 'dev'
                    # to check that it is a valid hash, make sure it is the right
                    # length and made up of valid hash characters
                    if re.match('[a-fA-F0-9]{40}$', sa['commit']):
                        return sa['commit']
        # again, default to setting this to release
        return 'dev' #'release'


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


