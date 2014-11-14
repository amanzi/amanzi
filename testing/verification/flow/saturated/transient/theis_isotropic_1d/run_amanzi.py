import os, sys, subprocess

def run_amanzi(input_file, directory=None):
    if directory is None:
        directory = os.getcwd()

    CWD =os.getcwd()
    run_directory= os.path.join(CWD,"amanzi-output")

    if os.path.isdir(run_directory):
        return
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
        if (os.environ['AMANZI_RUN_PARALLEL']):
            mpi_nprocs=os.environ['AMANZI_MPI_MAXPROCS']
            mpi_cmd = os.environ['AMANZI_MPI_EXEC'] + " " + os.environ['AMANZI_MPI_NUMPROCS_FLAG']+" "+mpi_nprocs+" "
            ierr = subprocess.call([mpi_cmd + executable + " --xml_file="+input_file + " --xml_schema="+xmlschema_path], stdout=stdout_file, stderr= subprocess.STDOUT, shell=True)
        else:
            ierr = subprocess.call([executable, "--xml_file="+input_file, "--xml_schema="+xmlschema_path], stdout=stdout_file, stderr= subprocess.STDOUT)
        #endif
        
    finally:
        os.chdir(CWD)

        return ierr
