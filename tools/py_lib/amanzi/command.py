################################################################################

import os, sys
import shlex

################################################################################

class CommandInterface:
    
    def __init__(self,command,args=None):
        
        self.command = command
        self.args = []
        self.exit_code=0
        self.output=''

        if args is not None:
            self._parse_arg_list(args)

        try:
            import subprocess
            self.use_ospipe = False
        except ImportError:
            self.use_ospipe = True
        except Exception:
            print("Unexpected error:",sys.exec_info()[0])
            raise
            

    def _dump_state(self):
        print('command=',self.command)
        print('args=',self.args)
        print('exit_code=',self.exit_code)
        print('output=',self.output)
        print('use_ospipe=',self.use_ospipe)

    def _parse_arg_list(self,args):
        list_args=[]
        if isinstance(args,list):
            list_args = shlex.split(" ".join(args))
        elif isinstance(args,str):
            list_args = shlex.split(args)
        else:
            raise TypeError('args must be of type list or str')
        self.args = list_args 
        return list_args

    def _build_run_command(self):
        return self.command + ' ' + ' '.join(self.args)

    def _parse_shell_exit(self,pattern):
        m = pattern.findall(self.output)
        if m is not None:
            idx = len(m) - 1
            self.exit_code = m[idx]

        return self.exit_code    

    def _remove_shell_exit(self,pattern):
        self.output = pattern.sub('',self.output)
       
    def set_args(self,args):
        self.args = self._parse_arg_list(args)
        return self.args

    def add_args(self,args):
        if len(self.args) == 0:
            self.set_args(args)
        else:  
            new_args = self._parse_arg_list(args)
            self.args.append(new_args)
        return self.args

    def search_args(self,target,index=None):
        index = -1
        if target in self.args:
            index = self.args.index(target)

        return index

    def clear_args(self):
        self.args = []
        return 

    def run(self):
        if self.use_ospipe:
            self._ospipe_run()
        else:
            self._subprocess_run()

        return self.exit_code    
    
    def _subprocess_run(self):

        import subprocess
        from subprocess import Popen,PIPE,STDOUT
        import time
        import signal

        run_command=self._build_run_command()
        try:
            pipe = Popen(run_command,shell=True,stdout=PIPE,stderr=STDOUT)
        except ValueError:
            raise ValueError('Incorrect arguments in Popen')
        except Exception:
            raise
        else:
            output=pipe.stdout
            self.ouput=output.read()
            output.close

        # Do not leave until the process completes!
        while pipe.poll() is None:
            pass

        # Set the return code
        self.exit_code=pipe.returncode 
        
        # Dump out the output
        print(self.output)

        return self.exit_code

    def _ospipe_run(self):
        import re
        run_command = self._build_run_command()

        # Need the '$' to delete the last print just
        # in case the command output also has this 
        # print out!
        pattern = re.compile('SHELL_EXIT=(\\d+)$')
        
        run_command = run_command + '; echo SHELL_EXIT=$?'
        (child_stdin, child_outerr) =os.popen4(run_command)
        child_stdin.close()
        self.output = child_outerr.read()
        self.exit_code = child_outerr.close()
        if self.exit_code is None:
            self._parse_shell_exit(pattern)
        self._remove_shell_exit(pattern)    
            

        return self.exit_code

################################################################################

def Command(command,args=None):
    cmd = CommandInterface(command,args)
    cmd.run()
    
    return cmd

################################################################################
if __name__ == '__main__':

    command='ls'
    ci=CommandInterface(command)
    ci.set_args('-lat')
    ci._dump_state()

    args=['-l', '-a']
    cmd = Command(command,args=args)
    cmd._dump_state()

    yacmd = Command('which','h5diff')

    bad_cmd=Command('dummy_exe')
    bad_cmd._dump_state()

