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

        if args != None:
            self._parse_arg_list(args)

        try:
            import subprocess
            self.use_ospipe = False
        except:
            self.use_ospipe = True

    def _dump_state(self):
        print 'command=',self.command
        print 'args=',self.args
        print 'exit_code=',self.exit_code
        print 'output=',self.output

    def _parse_arg_list(self,args):
        list_args=[]
        if isinstance(args,list):
            list_args = shlex.split(" ".join(args))
        elif isinstance(args,str):
            list_args = shlex.split(args)
        else:
            print 'Unknown instance that is not a string or list'
            sys.exit(1)

        return list_args    
       
    def set_args(self,args):
        self.args = self._parse_arg_list(args)
        return self.args

    def add_args(self,args):
        new_args = self._parse_arg_list(args)
        for item in new_args:
            self.args.append(item)

        return self.args

    def run(self):
        if self.use_ospipe == True:
            self._ospipe_run()
        else:
            self._subprocess_run()

        return self.exit_code    
    
    def _subprocess_run(self):

        try:
            import subprocess
            from subprocess import Popen,PIPE,STDOUT
            run_command=self.command + ' ' + ' '.join(self.args)
            pipe = Popen(run_command,shell=True,stdout=PIPE,stderr=STDOUT)
            output = pipe.stdout
            self.output = output.read()
            output.close()
            self.exit_code = pipe.wait()
            print self.output
        except:
             print 'Python module subprocess is not available'
             sys.exit(1)

        return self.exit_code

    def _ospipe_run(self):
        print 'Not support at this time'
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



    bad_cmd=Command('dummy_exe')
    bad_cmd._dump_state()
