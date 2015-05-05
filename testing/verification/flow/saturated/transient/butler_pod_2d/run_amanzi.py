import os, sys, subprocess

def run_amanzi(input_file, directory=None):
    if directory is None:
        directory = os.getcwd()

    CWD =os.getcwd()
    run_directory= os.path.join(CWD,"amanzi-output")
    run_stdout=os.path.join(run_directory,"stdout.out")

    if os.path.isdir(run_directory):
        if os.path.isfile(run_stdout):
            if ( "Amanzi::SIMULATION_SUCCESSFUL" in open(run_stdout).read() ):
                return
    else:
        os.mkdir(run_directory) 

    os.chdir(run_directory)

    # ensure any required files are present
    os.symlink("../mesh_cylinder.exo","mesh_cylinder.exo")
    
    # ensure that Amanzi's executable exists
    try:
        executable = os.environ['AMANZI_EXE']
    except KeyError:
        raise RuntimeError("Missing Amanzi installation, please set the AMANZI_EXE environmental variable.")

    if not os.path.isfile(executable):
        raise RuntimeError("The executable,\n"+executable+"\nis missing. Please build and install Amanzi.")

    # ensure that Amanzi's XML schema file exists
    try:
        xmlschema_file=os.environ['AMANZI_XSD']
    except KeyError:
        raise RuntimeError("Missing XML schema file, please set the AMANZI_XSD environment variable.")

    # determine if this is a serial or parallel run
    try:
        run_parallel_env=os.environ['AMANZI_RUN_PARALLEL'].lower()
        if (run_parallel_env=="yes" or run_parallel_env=="true" or run_parallel_env=="1"):
            run_parallel=True
        else:
            run_parallel=False
    except KeyError:
        run_parallel=False

    try:
        stdout_file = open("stdout.out", "w")
        if (run_parallel):
            mpi_nprocs=os.environ['AMANZI_MPI_MAXPROCS']
            mpi_cmd = os.environ['AMANZI_MPI_EXEC'] + " " + os.environ['AMANZI_MPI_NUMPROCS_FLAG']+" "+mpi_nprocs+" "
            mpi_job = mpi_cmd + executable + " --xml_file=" + input_file + " --xml_schema=" + xmlschema_file
            #stdout_file.write("Executing: "+ mpi_job + "\n" )
            ierr = subprocess.call([mpi_job], stdout=stdout_file, stderr= subprocess.STDOUT, shell=True)
        else:
            ser_job = [executable, "--xml_file="+input_file, "--xml_schema="+xmlschema_file]
            #print "Executing: ", ser_job
            ierr = subprocess.call(ser_job, stdout=stdout_file, stderr= subprocess.STDOUT)
        #endif
        
    finally:
        os.chdir(CWD)

        return ierr
