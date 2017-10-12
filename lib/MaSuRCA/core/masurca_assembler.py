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

from MaSuRCA.core.masurca_utils import masurca_utils

from DataFileUtil.DataFileUtilClient import DataFileUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
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


class MaSuRCA_Assembler(object):
    INVALID_WS_OBJ_NAME_RE = re.compile('[^\\w\\|._-]')
    INVALID_WS_NAME_RE = re.compile('[^\\w:._-]')

    PARAM_IN_CS_NAME = 'output_contigset_name'
    PARAM_IN_READ_LIB = 'read_libraries'
    PARAM_IN_JUMP_LIB = 'jump_libraries'
    PARAM_IN_PACBIOS = 'pacbios'
    PARAM_IN_OTHER = 'read_libraries'
    MaSuRCAR_PROJECT_DIR = 'masurca_project_dir'
    MaSuRCA_OUT_DIR = 'MaSuRCA_Output'

    def __init__(self, config, provenance):
        """
        __init__
        """
        # BEGIN_CONSTRUCTOR
        self.workspace_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.provenance = provenance

        self.au = AssemblyUtil(self.callback_url)

        self.scratch = os.path.join(config['scratch'], str(uuid.uuid4()))
        _mkdir_p(self.scratch)

        self.masurca_version = 'MaSuRCA-' + os.environ['M_VERSION']
        self.proj_dir = self.create_proj_dir(self.scratch)
        self.m_utils = masurca_utils(self.proj_dir, config)

        # from the provenance, extract out the version to run by exact hash if possible
        self.my_version = 'release'
        if len(provenance) > 0:
            if 'subactions' in provenance[0]:
                self.my_version = self.get_version_from_subactions('kb_MaSuRCA', provenance[0]['subactions'])
        print('Running kb_MaSuRCA version = ' + self.my_version)
        # END_CONSTRUCTOR
        pass


    def run_masurca_app2(self, params):
        # 1. validate & process the input parameters
        validated_params = self.m_utils.validate_params(params)

        wsname = params['workspace_name']

        # 2. create the configuration file 
        config_file = self.m_utils.construct_masurca_assembler_cfg(validated_params)

        # 3. run masurca against the configuration file to generate the assemble.sh script 
        if os.path.isfile(config_file):
            assemble_file = self.m_utils.generate_assemble_script(config_file)

        # 4. run the assemble.sh script to do the heavy-lifting
        if os.path.isfile(assemble_file):
            assemble_ok = self.m_utils.run_assemble(assemble_file)
        else:
            assemble_ok = -1

        # 5. report the final results
        returnVal = {
            "report_ref": None,
            "report_name": None
        }

        # 6. save the assembly to KBase and, if everything has gone well, create a report
        ca_dir = os.path.join(self.proj_dir, 'CA*')
        contig_fa_file = 'dedup.genome.scf.fasta'
        contig_fa_file = os.path.join(ca_dir, contig_fa_file)

        #if (assemble_ok == 0 and os.path.isfile(contig_fa_file)):
        fa_file_path = self.find_file_path(self.proj_dir, contig_fa_file)
        log("Found contig_fa_file with path {}".format(fa_file_path))
        if (assemble_ok == 0 and fa_file_path != ''):
            fa_file_path = os.path.join(self.proj_dir, fa_file_path)
            self.m_utils.save_assembly(fa_file_path, wsname, params[self.PARAM_IN_CS_NAME])
            if params['create_report'] == 1:
                report_name, report_ref = self.m_utils.generate_report(contig_fa_file, params, ca_dir, wsname)
                returnVal = {'report_name': report_name, 'report_ref': report_ref}
        else:
            log("run_assemble process failed.")

        return returnVal

    def find_file_path(self, search_dir, search_file_name):
        for dirName, subdirList, fileList in os.walk(search_dir):
            log('Found directory: {}'.format(dirName))
            for fname in fileList:
                log('Found file: {}'.format(fname))
                if fname == search_file_name:
                    return os.path.join(dirName, fname)
        return ''


    def run_masurca_app(self, params):
        # 1. validate & process the input parameters
        validated_params = self.m_utils.validate_params(params)

        wsname = params['workspace_name']

        # 2. create the configuration file 
        config_file = self.m_utils.construct_masurca_config(validated_params)

        # 3. run masurca against the configuration file to generate the assemble.sh script 
        if os.path.isfile(config_file):
            assemble_file = self.m_utils.generate_assemble_script(config_file)

        # 4. run the assemble.sh script to do the heavy-lifting
        if os.path.isfile(assemble_file):
            assemble_ok = self.m_utils.run_assemble(assemble_file)
        else:
            assemble_ok = -1

        # 5. report the final results
        returnVal = {
            "report_ref": None,
            "report_name": None
        }

        # 6. save the assembly to KBase and, if everything has gone well, create a report
        ca_dir = os.path.join(self.proj_dir, 'CA')
        contig_fa_file = 'dedup.genome.scf.fasta'
        contig_fa_file = os.path.join(ca_dir, contig_fa_file)

        if (assemble_ok == 0 and os.path.isfile(contig_fa_file)):
            self.m_utils.save_assembly(contig_fa_file, wsname, params[self.PARAM_IN_CS_NAME])
            if params['create_report'] == 1:
                report_name, report_ref = self.m_utils.generate_report(contig_fa_file, params, ca_dir, wsname)
                returnVal = {'report_name': report_name, 'report_ref': report_ref}
        else:
            log("run_assemble process failed.")

        return returnVal


    def create_proj_dir(self, home_dir):
        """
        creating the project directory for MaSuRCA
        """
        prjdir = os.path.join(home_dir, self.MaSuRCAR_PROJECT_DIR)
        _mkdir_p(prjdir)
        self.proj_dir = prjdir

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

