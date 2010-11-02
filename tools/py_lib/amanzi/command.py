################################################################################

import os, sys
import shlex

################################################################################

class CommandInterface:
    
    def __init__(self,command,args=None):
        
        self.command = command
        if args != None:
            if isinstance(args,list):
                self.args=shlex.split(" ".join(args))
            elif isinstance(args,str):
                self.args=shlex.split(args)
            else:
                print 'Invalid argument type'
                sys.exit(1)
        else:
            self.args=[]

        self.exit_code=0
        self.output=''

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

    def subprocess_run(self):

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
            print 'RETURN CODE',pipe.returncode
        except:
             print 'Python module subprocess is not available'
             sys.exit(1)

        return self.exit_code

    def ospipe_run(self):
        print 'Not support at this time'
        return self.exit_code

def Command(command,args=None):
    cmd = CommandInterface(command,args)
    if cmd.use_ospipe != True:
        cmd.subprocess_run()
    else:
        cmd.ospipe_run()
    return cmd



if __name__ == '__main__':

    command='ls'
    args='-last'
    ci=CommandInterface(command,args)
    ci._dump_state()

    args=['-l', '-a']
    cmd = Command(command,args=args)
    cmd._dump_state()

    bad_cmd=Command('dummy_exe')
    bad_cmd._dump_state()
