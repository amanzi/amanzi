######### Running Amanzi to generate an output file #############
#
# This script will create a run directory, copy specified files
# to it and run the simulation from it. The log is saved in file
# stdout.out. Amanzi-S runs from the base directory.
#
# input_file = XML input file located in the base directory.
#              It is not copied to the run directory.
#
# man_np = The maximum number of cores to use for parallel run.
#          It must be specified via the enviromental variable, e.g.
#          export AMANZI_MPI_NP=10
#
# copyfiles = The list of files to copy from the base directory
#             to the run directory. The file names should not
#             contain path in them, e.g. 
#             copyfile=[mesh.exo, tritium.gdb].
# 
# subdirectory = The name for the run directory. Default name is 
#                amanzi-output. This name should not have path
#                in it.
#
#################################################################
import os, sys, subprocess
import shutil

def run_amanzi(input_file, max_np, copyfiles=None, subdirectory=None):
    cwd = os.getcwd()

    # set up the run directory
    if subdirectory is None:
       run_directory = os.path.join(cwd, "output")
    else:
       run_directory = os.path.join(cwd, subdirectory)

    # set up the log file and verify it content
    run_stdout = os.path.join(run_directory, "stdout.out")

    print("  trying ", input_file)
    if os.path.isdir(run_directory):
       if os.path.isfile(run_stdout):
          if ( "Amanzi::SIMULATION_SUCCESSFUL" in open(run_stdout).read() ):
             print("    data exist, skipping the test") 
             os.chdir(cwd)
             return
    else:
        os.mkdir(run_directory) 

    # copy files to the run directory
    os.chdir(run_directory)

    if copyfiles is None:
       pass
    else:
       for cf in copyfiles:
           oldfile = "../" + cf
           newfile = run_directory + "/" + cf
           shutil.copyfile(oldfile, newfile)

    # miscalleneous
    xml_cmd = "--xml_file="

    # ensure that Amanzi's executable exists
    try:
        path = os.path.join(os.environ['AMANZI_INSTALL_DIR'], 'bin')
    except KeyError:
        raise RuntimeError("Missing Amanzi installation, check environmental variable AMANZI_INSTALL_DIR.")
    executable = os.path.join(path, "amanzi")

    if not os.path.isfile(executable):
       raise RuntimeError("Missing Amanzi installation, please build and install Amanzi.")

    # run the simulator
    try:
        logfile = open("stdout.out", "w")
        mpi_exec = os.getenv('AMANZI_MPI_EXEC', 'mpirun')
        mpi_np = os.getenv('AMANZI_MPI_NP', '1')
        mpi_np = str(min(int(mpi_np), max_np))
        if (mpi_np != '1'):
           mpi_cmd = mpi_exec + " -np " + mpi_np + " --oversubscribe " + executable + " " + xml_cmd + input_file
           ierr = subprocess.call(mpi_cmd, stdout=logfile, stderr=subprocess.STDOUT, shell=True)
        else:
           args = xml_cmd + input_file
           ierr = subprocess.call([executable, args], stdout=logfile, stderr=subprocess.STDOUT)

    finally:
        pass
        
    os.chdir(cwd)

    return ierr
