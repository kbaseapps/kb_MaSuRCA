import re
import time
import json
import os
from pprint import pprint, pformat

from MaSuRCA.core.Program_Runner import Program_Runner
from DataFileUtil.DataFileUtilClient import DataFileUtil
from Workspace.WorkspaceClient import Workspace as Workspace
from KBaseReport.KBaseReportClient import KBaseReport
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil

from file_util import (
    valid_string,
    check_reference,
    get_unique_names,
    fetch_fasta_from_object,
    fetch_reads_refs_from_sampleset,
    fetch_reads_from_reference,
    extract_geneCount_matrix
)

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))

class masurca_utils:

    def __init__(self, scratch_dir, workspace_url, callback_url, srv_wiz_url, provenance):
        self.workspace_url = workspace_url
        self.callback_url = callback_url
        self.srv_wiz_url = srv_wiz_url
        self.au = AssemblyUtil(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url, service_ver='beta')
        self.scratch = scratch_dir
        self.working_dir = scratch_dir
        self.prog_runner = Program_Runner(self.STAR_BIN, self.scratch)
        self.provenance = provenance
        self.ws_client = Workspace(self.workspace_url)

        self.parallel_runner = KBParallel(self.callback_url)
        self.qualimap = kb_QualiMap(self.callback_url, service_ver='dev')
        self.set_api_client = SetAPI(self.srv_wiz_url, service_ver='dev')
        self.eu = ExpressionUtils(self.callback_url, service_ver='beta')

    def _replaceSectionText(self, orig_txt, begin-patn, end-patn, repl_txt, filename):
        """
        replace a section of text between lines begin-patn and end-patn with repl_text in file named by filename
        """
        ret_val = false
        try:
            # create regular expression pattern
            repl = re.compile('begin-patn.*?#end-patn', re.DOTALL)

            # open file
            f = open(filename, 'r')
            txt = f.read()
            f.close()

            # chop text between #chop-begin and #chop-end
            txt_replaced = repl.sub(repl_txt, txt)
            pprint(txt_replaced)

            # save result
            f = open(filename, 'w')
            f.write(txt_replaced)
            f.close()
        except ValueError as ve:
            log('File modification raised error:\n')
            pprint(ve)
        else: #no exception raised
            ret_val = true

        return ret_val


