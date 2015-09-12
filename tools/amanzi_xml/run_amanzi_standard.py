#### Running Amanzi to generate an output file ####
import os, sys, subprocess

def run_amanzi(input_file, max_np, directory=None):
    if directory is None:
        directory = os.getcwd()
    cwd = os.getcwd()
    run_directory = os.path.join(cwd, "amanzi-output")
    run_stdout = os.path.join(run_directory, "stdout.out")

    if os.path.isdir(run_directory):
        if os.path.isfile(run_stdout):
            if ( "Amanzi::SIMULATION_SUCCESSFUL" in open(run_stdout).read() ):
                print("  data exist, skipping the test") 
                return
    else:
        os.mkdir(run_directory) 

    os.chdir(run_directory)
    
    # ensure that Amanzi's executable exists
    try:
        path = os.path.join(os.environ['AMANZI_INSTALL_DIR'], 'bin')
    except KeyError:
        raise RuntimeError("Missing Amanzi installation, check environmental variable AMANZI_INSTALL_DIR.")
    executable = os.path.join(path, "amanzi")

    if not os.path.isfile(executable):
        raise RuntimeError("Missing Amanzi installation, please build and install Amanzi.")

    try:
        logfile = open("stdout.out", "w")
        mpi_exec = os.getenv('AMANZI_MPI_EXEC', 'mpirun')
        mpi_np = os.getenv('AMANZI_MPI_NP', '1')
        if (mpi_np != '1'):
           mpi_cmd = mpi_exec + " -np " + mpi_np + " " + executable + " --xml_file=../" + input_file
           ierr = subprocess.call(mpi_cmd, stdout=logfile, stderr=subprocess.STDOUT, shell=True)
        else:
           ierr = subprocess.call([executable, "--xml_file=../" + input_file], stdout=logfile, stderr=subprocess.STDOUT)
        
    finally:
        os.chdir(cwd)

    return ierr
