#### Running Amanzi to generate an output file ####
import os

def run_amanzi(input_file, directory=None):
    if directory is None:
        directory = os.getcwd()

    # ensure that Amanzi's executable exists
    try:
        path = os.path.join(os.environ['AMANZI_INSTALL_DIR'],'bin')
    except KeyError:
        raise RunTimeError("Missing Amanzi installation, please set the AMANZI_INSTALL_DIR environmental variable.")
    executable = os.path.join(path, "amanzi")
    if not os.path.isfile(executable):
        raise RunTimeError("Missing Amanzi installation, please build and install Amanzi.")

    # set up the run directory and cd into it
    run_directory = os.path.join(directory,"output")
    if os.path.isdir(run_directory):
        [os.remove(f) for f in os.listdir(run_directory)]
    else:
        os.mkdir(run_directory)
    os.chdir(run_directory)

    # run amanzi
    try:
        stdout = open(input_file.split(".")[0]+".log",'w')
        stderr = open(input_file.split(".")[0]+".err",'w')
        try:
            ierr = os.spawnl(os.P_WAIT, executable, ["--xml_file=%s"%input_file,], stdout=stdout)
        finally:
            stdout.close()
            stderr.close()
    finally:
        os.chdir(directory)

    return ierr
