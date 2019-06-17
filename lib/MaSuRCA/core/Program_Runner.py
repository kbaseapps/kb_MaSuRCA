import sys
import subprocess


class Program_Runner:

    def __init__(self, cmd, work_dir):
        self.work_dir = work_dir
        self.executableName = cmd

    def run(self, command, cwd_dir=None):
        '''
        options is an array of command-line parameters passed
        to the RQCFilter App
        '''
        cmmd = command
        if not cwd_dir:
            cwd_dir = self.work_dir

        res = subprocess.Popen(cmmd, cwd=cwd_dir, shell=False)
        exitCode = res.wait()

        if (exitCode == 0):
            print('\n', ' '.join(cmmd),
                  ' was executed successfully, exit code: ' + str(exitCode))
        else:
            if sys.exc_info()[0]:
                print('Error > ', sys.exc_info()[0])
            raise ValueError('Error running command: ' + ' '.join(cmmd) +
                             '\nExit Code: ' + str(exitCode))
        return exitCode
