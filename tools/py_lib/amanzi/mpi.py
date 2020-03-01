################################################################################

import sys, os
import types
from amanzi.command import CommandInterface

################################################################################

class MpiInterface(CommandInterface):

    def __init__(self,mpirun_exe='mpirun',args=None):
        
        # Check the execute command
        CommandInterface.__init__(self,mpirun_exe,args)

        self.mpirun_exe = mpirun_exe

          
    def num_procs(self,n):
        ntype = type(n)
        if isinstance(ntype,int):
            np = n
        elif isinstance(ntype,str):
            try:
                stripped = str(int(n))
            except Exception:
                print(n, 'is not an integer')

            np = int(n)

        self.np = np

        # Search the args to see if the number of procs has been set
        possible_np_args = [ '-n', '--n', '-np', '--np' ]
        n_try = 0
        max_try = len(possible_np_args)
        arg_index = -1
        while n_try < max_try and (arg_index < 0 ):
            opt = possible_np_args[n_try]
            arg_index =  self.search_args(opt)
            n_try = n_try + 1

        if arg_index >= 0:
            self.args[arg_index+1] = str(np)
        else:
            self.args.insert(0,str(np))
            self.args.insert(0,'-np')


    def _dump_state(self):
        print('')
        print('################################################################################')
        print('')
        print('command:', self.command)
        print('args:', self.args)
        print('exit_code:', self.exit_code)
        print('')
        print('################################################################################')
        print('')

    def run(self,binary=None,binary_args=None):
        if binary is not None:
            self.add_args(binary)

        if binary_args is not None:
            self.add_args(binary_args)

        CommandInterface.run(self)

        return self.exit_code
        
################################################################################
if __name__ == '__main__':

    mpi = MpiInterface()

    print(mpi.command)
    print(mpi.mpirun_exe)
    print(mpi.args)

    # Passing args as a list
    mpi.num_procs(4)
    mpi.run('hello_world',['-a', '--solver=jack'])
    mpi._dump_state()

    # Passing args as a string
    mpi.clear_args()
    mpi.num_procs(4)
    mpi.run('new_binary', '-a --preifx')
    mpi._dump_state()

    # This test should fail
    try:
        mpi.num_procs('blah')
        mpi._dump_state()
    except Exception:
        print('Passed the invalid proc test')

    # Resetting the number of procs
    mpi.num_procs(8)
    mpi._dump_state()

