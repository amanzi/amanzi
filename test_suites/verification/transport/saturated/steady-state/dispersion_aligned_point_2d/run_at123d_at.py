#### Running AT123D-AT to generate analytic solutions ####

import os, sys, shutil, subprocess, re

def run_at123d(input_file, directory=None, overwrite=False):

    if directory is None:
        directory = "at123d-at"

    CWD =os.getcwd()
    run_directory= os.path.join(CWD,directory)

    if not os.path.isdir(run_directory):
        os.mkdir(run_directory) 

    # Change to run directory
    os.chdir(run_directory)

    # Set the name of the log file
    at123d_setup=os.path.splitext(input_file)[0]+"_setup.out"
    at123d_soln=os.path.splitext(input_file)[0]+"_soln.out"

    # Return if the log file is already there.
    if ( os.path.isfile(at123d_setup) and os.path.isfile(at123d_soln) and not overwrite ):
        os.chdir(CWD)
        return

    # ensure that AT123D-AT executable exists
    if ( os.environ.get('AT123DAT_EXE') ):
        try:
            os.path.isfile(os.environ.get('AT123DAT_EXE'))
        except KeyError:
            raise RuntimeError("The executable specifed AT123DAT_EXE "+os.environ.get('AT123DAT_EXE')+" does not exist!")
    else:
       raise Exception("Please set the AT123DAT_EXE environmental variable to specify the full path and name of the executable")

    executable = os.environ.get('AT123DAT_EXE') 

    try:
        at123d_log=os.path.splitext(input_file)[0]+".log"
        stdout_file = open(at123d_log, "w")
        stdin_file = open(input_file,"r")
        ierr = subprocess.call([executable, ], stdin=stdin_file, stdout=stdout_file, stderr= subprocess.STDOUT)
        
    finally:
        os.chdir(CWD)
        return ierr
