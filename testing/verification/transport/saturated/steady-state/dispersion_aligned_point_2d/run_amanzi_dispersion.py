#### Running Amanzi to generate an output file ####
import os, sys, shutil, subprocess, re

def run_amanzi(input_file, directory=None, parameters=None, mesh_file=None, overwrite=False):

    if directory is None:
        directory = "amanzi-output"

    CWD =os.getcwd()
    run_directory= os.path.join(CWD,directory)

    if not os.path.isdir(run_directory):
        os.mkdir(run_directory) 

    # Change to run directory
    os.chdir(run_directory)

    # Obersvation file - should really come from the input file.
    obs_file="obs5_2_1_r5.out"

    # Return if there is already output and we aren't overwriting
    if ( os.path.isfile(input_file) and os.path.isfile(obs_file) and not overwrite ):
        os.chdir(CWD)
        return

    # copy file
    if parameters is None and mesh_file is None:
        shutil.copy(os.path.join("..",input_file), run_directory)
        parameters={}
    else:
        orig_file=open(os.path.join("..",input_file),'r')
        new_file=open(input_file,"w")
        InMesh=False
        for line in orig_file:
            if ( not InMesh ):
                if ( "Read Mesh File" in line ):
                    InMesh=True
                    new_file.write(line)
                else:
                    for name in parameters:
                        if ( name in line ): 
                            new_file.write(re.sub(r'value=.*"','value="'+parameters[name]+'"', line))
                        else:
                            new_file.write(line)
            elif ( "File" in line ):
                new_file.write(re.sub(r'value=.*"','value="'+mesh_file+'"', line))
                InMesh = False
            else:
                new_file.write(line)

        orig_file.close()
        new_file.close()


    # ensure that Amanzi's executable exists
    try:
        path = os.path.join(os.environ['AMANZI_INSTALL_DIR'],'bin')
    except KeyError:
        raise RuntimeError("Missing Amanzi installation, please set the AMANZI_INSTALL_DIR environmental variable.")
    executable = os.path.join(path, "amanzi")
    print executable

    if not os.path.isfile(executable):
        raise RuntimeError("Missing Amanzi installation, please build and install Amanzi.")

    try:
        stdout_file = open("stdout.out", "w")
        ierr = subprocess.call([executable, "--xml_file="+input_file], stdout=stdout_file, stderr= subprocess.STDOUT)
        
    finally:
        os.chdir(CWD)

        return ierr
