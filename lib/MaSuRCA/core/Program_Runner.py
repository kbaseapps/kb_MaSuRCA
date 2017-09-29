import os
import subprocess


class Program_Runner:

    def __init__(self, cmd, work_dir):
        self.work_dir = work_dir
        self.executableName = cmd

    def run(self, command, cwd_dir=None):
        ''' options is an array of command-line parameters passed to the RQCFilter App '''
        cmmd = command

        if not cwd_dir:
          cwd_dir = self.work_dir

        print('\nRunning: ' + ' '.join(cmmd))
        res = subprocess.Popen(cmmd, cwd=cwd_dir, shell=False)

        exitCode = res.wait()

        if (exitCode == 0):
            print('\n' + ' '.join(cmmd) + ' was executed successfully, exit code was: ' + str(exitCode))
        else:
            print "Error > ",sys.exc_info()[0]
            raise ValueError('Error running command: ' + ' '.join(cmmd) + '\n' +
                             'Exit Code: ' + str(exitCode))

        return exitCode

    def runi_1(self, command, cwd_dir=None):
        ''' options is an array of command-line parameters passed to the RQCFilter App '''
        cmmd = command

        if not cwd_dir:
          cwd_dir = self.work_dir

        print('\nRunning: ' + ' '.join(cmmd))
        exitCode = -1
        try:
            res = subprocess.Popen(cmmd, cwd=cwd_dir, shell=False)
            output,error = res.communicate()
            exitCode = res.returncode
            if output:
                print "ret> ",res.returncode
                print "OK> output ",output
            if error:
                print "ret> ",res.returncode
                print "Error> error ",error.strip()
        except OSError as e:
            print "OSError > ",e.errno
            print "OSError > ",e.strerror
            print "OSError > ",e.filename
        except:
            print "Error > ",sys.exc_info()[0]

        return exitCode

