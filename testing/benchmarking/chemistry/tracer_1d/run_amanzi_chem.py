#### Running Amanzi to generate an output file ####
#### Takes in chemistry database for native Amanzi chemistry
import os, sys, subprocess, shutil

def run_amanzi_chem(xml_input_file, run_path=None, chemfiles=None, directory=None):

    # get current directory
    if directory is None:
        directory = os.getcwd()
    CWD =os.getcwd()

    # run directory: provided or guessed
    if run_path is None:
        run_directory= os.path.join(CWD,"amanzi-output")
    else:
        run_directory= os.path.join(CWD,run_path)

    # is run directory populated already - do not rerun if so
    if os.path.isdir(run_directory):
        if os.listdir(run_directory):
           return
    else:
        os.mkdir(run_directory) 

    # copy bgd file for native chemistry: needs to be in run directory
    if chemfiles is None:
        pass
    else:
        for chemfile in chemfiles:
           new_chemfile="{0}/{1}".format(run_directory, chemfile)
           shutil.copyfile(chemfile, new_chemfile)

    # change into run directory
    if run_path[0:7] == 'struct_':
       pass
    else:  
       os.chdir(run_directory)
    
    # ensure that Amanzi's executable exists
    try:
        path = os.path.join(os.environ['AMANZI_INSTALL_DIR'],'bin')
    except KeyError:
        raise RuntimeError("Missing Amanzi installation, please set the AMANZI_INSTALL_DIR environmental variable.")
    executable = os.path.join(path, "amanzi")

    if not os.path.isfile(executable):
        raise RuntimeError("Missing Amanzi installation, please build and install Amanzi.")

    try:
        stdout_file = open("stdout.out", "w")
        print("TEST: tracer")
        ierr = subprocess.call([executable, "--xml_file="+xml_input_file], stdout=stdout_file, stderr= subprocess.STDOUT)

    finally:
        os.chdir(CWD)

        return ierr
