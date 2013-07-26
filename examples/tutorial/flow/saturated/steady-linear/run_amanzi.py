#### Running Amanzi to generate an output file ####
import os, sys, subprocess

def run_amanzi(input_file, directory=None):
    if directory is None:
        directory = os.getcwd()
    CWD =os.getcwd()
    run_directory= os.path.join(CWD,"amanzi-output")

    if os.path.isdir(run_directory):
        answer = raw_input("Do you wish to run amanzi and delete the current observation file? (y or n): ")
        if answer == 'y':
            [os.remove(os.path.join(run_directory,f)) for f in os.listdir(run_directory)]
        elif answer =='n':
            return "Amanzi will not run!!"
    else:
        os.mkdir(run_directory) 
    os.chdir(run_directory)
    
    # ensure that Amanzi's executable exists
    try:
        path = os.path.join(os.environ['AMANZI_INSTALL_DIR'],'bin')
    except KeyError:
        raise RunTimeError("Missing Amanzi installation, please set the AMANZI_INSTALL_DIR environmental variable.")

    executable = os.path.join(path, "amanzi")

    if not os.path.isfile(executable):
        raise RunTimeError("Missing Amanzi installation, please build and install Amanzi.")

    try:
        stdout_file = open("stdout.out", "w")
        ierr = subprocess.call([executable, "--xml_file="+input_file], stdout=stdout_file, stderr= subprocess.STDOUT)
        
    finally:
        os.chdir(CWD)

    return ierr
