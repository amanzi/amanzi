#!/bin/sh
################################################################################
#                                                                              #
# Amanzi Build Script                                                          #
#                                                                              #
################################################################################

# ---------------------------------------------------------------------------- #
# Initialize                  
# ---------------------------------------------------------------------------- #

# Useful environment variabes
export PROJECT_DIR=/project/projectdirs
export ASCEM_GROUP=m1012
export ASCEM_HOME=${PROJECT_DIR}/${ASCEM_GROUP}

# Variables to improve code readability
TRUE=1
FALSE=0

# Initialize the module command
if [[ $MODULESHOME ]]; then
    . $MODULESHOME/init/bash
fi

# Amanzi source directory
amanzi_source_dir=`pwd`

# Number of parallel make threads
make_np=1

# Install directory
amanzi_install_dir=$HOME/amanzi

# Default programing environment
prg_env=gnu

# Default TPL CMake script file
init_file=

# Enable structured mesh
enable_structured=FALSE

# Enable unstructured mesh
enable_unstructured=FALSE

# Clobber the CMake build files before running cmake
clobber=FALSE

# CMake build type
amanzi_build_type=Debug

# Amanzi branch type
amanzi_branch=default

# Amanzi production install
amanzi_production=FALSE

# Amanzi enable HYPRE interface
amanzi_enable_hypre=FALSE

# Build as a batch job
batch_build=FALSE
batch_build_file=`pwd`/amanzi-build.pbs
batch_build_args=

# ---------------------------------------------------------------------------- #
# ASCEM Environment Variables                  
# ---------------------------------------------------------------------------- #
export ASCEM_HOME=/project/projectdirs/m1012
export ASCEM_TPL_INSTALL_DIR=${ASCEM_HOME}/tpls/install/${NERSC_HOST}
export ASCEM_MODULESHOME=${ASCEM_HOME}/modulefiles/${NERSC_HOST}

# ---------------------------------------------------------------------------- #
# Useful functions                  
# ---------------------------------------------------------------------------- #
exit_now()
{
    echo $1
    exit $2
}
amanzi_build_usage()
{
echo '
Usage: '"$0"' [options] DIR
where Amanzi source is located in DIR [default: '"${amanzi_source_dir}"']
Options: [Defaults in brackets]

    --help                  print this message
    --parallel=n            number of parallel make threads ['"${make_np}"']
    --install=DIR           install Amanzi in DIR ['"${amanzi_install_dir}"']
                            If DIR='production' then Amanzi will be installed
	                    in the ASCEM install location. Script will also run
			    the test suite and will not install if the test suite
			    fails. 
    --enable-structured     enable structured mesh code ['"${enable_structured}"']
    --enable-unstructured   enable structured mesh code ['"${enable_unstructured}"']
    --enable-hypre          enable HYPRE interface. Trilinos install MUST have HYPRE enabled ['"${amanzi_enable_hypre}"']
    --build-type=STRING     amanzi build type where STRING=Debug or Release ['"${amanzi_build_type}"']
    --init-file=FILE        initialize the CMake cache with FILE
    --prg-env=STRING        load PrgEnv module STRING ['"${prg_env}"']         
    --branch=STRING         build Amanzi branch STRING ['"${amanzi_branch}"']
    --clobber               remove CMakeCache.txt file and clean the build directory before running cmake ['"${clobber}"']
    --batch                 submit this build script as a batch job
'
  exit_now
}
print_info()
{
echo '
--------------------------------------------------------------------------------
amanzi_source_dir   ='"${amanzi_source_dir}"'
amanzi_install_dir  ='"${amanzi_install_dir}"' (production install: '"${amanzi_production}"')
amanzi_build_type   ='"${amanzi_build_type}"'
amanzi_branch       ='"${amanzi_branch}"'

prg_env             ='"${prg_env}"'
init_file           ='"${init_file}"'
clobber             ='"${clobber}"'
make_np             ='"${make_np}"'

enable_structured   ='"${enable_structured}"'
enable_unstructured ='"${enable_unstructured}"'
enable_hypre        ='"${amanzi_enable_hypre}"'

branch_build        ='"${batch_build}"'
branch_build_file   ='"${batch_build_file}"'
branch_build_args   ='"${branch_build_args}"'
--------------------------------------------------------------------------------
'
}
fix_build_type()
{
    dummy=`echo $1 | tr "[A-Z]" "[a-z]"`
    case $dummy in
	debug)
	    echo "Debug"
	    ;;
	release)
	    echo "Release"
	    ;;
	*)
	    echo "Error unknown build type $1"
	    exit 20
	    ;;
    esac
}
pe_env_tolower()
{
    echo ${PE_ENV} | tr "[A-Z]" "[a-z]"
}
get_compiler_version()
{
    version_var=${PE_ENV}_VERSION
    echo ${!version_var}
}
get_mpi_version()
{
    echo ${CRAY_MPICH2_VERSION}
}
mercurial_branch()
{
    echo `hg id --branch $1`
}
mercurial_global_id()
{
    echo `hg id --id $1`
}
mercurial_local_id()
{
    echo `hg id --num $1`
}
write_dflt_batch_file()
{
    root_name=amanzi-build-${prg_env}-${amanzi_build_type}
cat > $1  <<EOF
#PBS -S /bin/bash
#PBS -q regular
#PBS -l walltime=02:00:00
#PBS -N ${root_name}
#PBS -e ${root_name}.\$PBS_JOBID.err
#PBS -o ${root_name}.\$PBS_JOBID.out
#PBS -V

cd \$PBS_O_WORKDIR
$0 ${batch_build_args} ${amanzi_source_dir}

exit
EOF
}
submit_batch_script()
{
    echo "Submitting $1 as a batch script"
    qsub $1
    if [ $? -ne 0 ]; then
	exit_now "Failed to submit batch script" 50
    fi
}
production_install_dir()
{
    mpi_version="`get_mpi_version`"
    compiler_version="`get_compiler_version`"
    compiler_name="`pe_env_tolower`"
    #echo "${ASCEM_HOME}/install/${compiler_name}/mpich-${mpi_version}-${compiler_name}-${compiler_version}"
    echo "${HOME}/install-test/${compiler_name}/mpich-${mpi_version}-${compiler_name}-${compiler_version}"
} 
production_binary_name()
{
    source_dir=$1
    id="`mercurial_global_id $source_dir`"
    branch="`mercurial_branch $source_dir`"
    echo "amanzi_${id}_${branch}"
}
#LPRITCHload_prg_env()
#LPRITCH{
#LPRITCH    compiler_loaded=false
#LPRITCH    if [[ $PE_ENV ]]; then
#LPRITCH	compiler_loaded=`echo $PE_ENV | tr "[A-Z]" "[a-z]"`
#LPRITCH    fi
#LPRITCH
#LPRITCH    if [[ $compiler_loaded ]]; then
#LPRITCH	echo "Unloading ${compiler_loaded}"
#LPRITCH	module unload PrgEnv-${compiler_loaded}
#LPRITCH    fi
#LPRITCH
#LPRITCH    echo "Loading PrgEnv-$1"
#LPRITCH    #module use ${ASCEM_MODULESHOME}
#LPRITCH    #module load AmanziEnv-$1
#LPRITCH    module load PrgEnv-$1
#LPRITCH}
default_init_file()
{
    mpi_version="`get_mpi_version`"
    compiler_version="`get_compiler_version`"
    compiler_name="`pe_env_tolower`"
    root_install_dir=$ASCEM_TPL_INSTALL_DIR/mpich-${mpi_version}-${compiler_name}-${compiler_version}/${amanzi_build_type}
    init_file=${root_install_dir}/share/cmake/amanzi-tpl-config.cmake
}
update_repository_branch()
{
    save_dir=`pwd`
    echo "Updating  $1 to branch $2"
    cd $1
    if [ $? -ne 0 ]; then
	exit_now "Failed to change to directory $1 to update the branch"
    fi
    hg update $2
    if [ $? -ne 0 ]; then
	exit_now "Failed to update $1 to branch $2"
    fi
    cd $save_dir
}
# ---------------------------------------------------------------------------- #
# Process Command Line Options                  
# ---------------------------------------------------------------------------- #
amanzi_source_dir=`pwd`
for a in "$@"; do
    a_is_valid=FALSE
    if echo $a | grep "^--clobber" > /dev/null 2> /dev/null; then
	clobber=TRUE
	a_is_valid=TRUE
	batch_build_args="${batch_build_args} $a"
    fi
    if echo $a | grep "^--batch" > /dev/null 2> /dev/null; then
	batch_build=TRUE
	a_is_valid=TRUE
	if [ ! -z "${PBS_JOBID}" ]; then
	    exit_now "--batch is not allowed in a batch job" 10
	fi
    fi
    if echo $a | grep "^--branch=" > /dev/null 2> /dev/null; then
	amanzi_branch=`echo $a | sed "s/^--branch=//"`
	a_is_valid=TRUE
	batch_build_args="${batch_build_args} $a"
    fi
    if echo $a | grep "^--install=" > /dev/null 2> /dev/null; then
	amanzi_install_dir=`echo $a | sed "s/^--install=//"`
	a_is_valid=TRUE
	batch_build_args="${batch_build_args} $a"
    fi
    if echo $a | grep "^--parallel=" > /dev/null 2> /dev/null; then
	make_np=`echo $a | sed "s/--parallel=//"`
	a_is_valid=TRUE
	batch_build_args="${batch_build_args} $a"
    fi
    if echo $a | grep "^--enable-structured" > /dev/null 2> /dev/null; then
	enable_structured=TRUE
	a_is_valid=TRUE
	batch_build_args="${batch_build_args} $a"
    fi
    if echo $a | grep "^--enable-unstructured" > /dev/null 2> /dev/null; then
	enable_unstructured=TRUE
	a_is_valid=TRUE
	batch_build_args="${batch_build_args} $a"
    fi
    if echo $a | grep "^--enable-hypre" > /dev/null 2> /dev/null; then
	amanzi_enable_hypre=TRUE
	a_is_valid=TRUE
	batch_build_args="${batch_build_args} $a"
    fi
    if echo $a | grep "^--init-file=" > /dev/null 2> /dev/null; then
	init_file=`echo $a | sed "s/^--init-file=//"`
	a_is_valid=TRUE
	batch_build_args="${batch_build_args} $a"
    fi
    if echo $a | grep "^--prg-env=" > /dev/null 2> /dev/null; then
	prg_env=`echo $a | sed "s/^--prg-env=//"`
	a_is_valid=TRUE
	batch_build_args="${batch_build_args} $a"
    fi
    if echo $a | grep "^--build-type=" > /dev/null 2> /dev/null; then
	build_type=`echo $a | sed "s/^--build-type=//"`
	amanzi_build_type="`fix_build_type $build_type`"
	a_is_valid=TRUE
	batch_build_args="${batch_build_args} $a"
    fi
    if echo $a | grep "^--help" > /dev/null 2> /dev/null; then
	amanzi_build_usage
    fi
    if [[ -d $a ]]; then
	if [[ -e "$a/.hg" ]]; then
	  amanzi_source_dir=$a
	  a_is_valid=TRUE
        else
	  echo "$a is not a valid mercurial repository"
	  a_is_valid=FALSE
        fi  
    fi
    if [ "${a_is_valid}" == "FALSE" ]; then
	echo "'$a' is an invalid option or an invalid source directory"
	amanzi_build_usage
    fi
done

echo "Building Amanzi source located in $amanzi_source_dir"

# ---------------------------------------------------------------------------- #
# Lock file control                  
# ---------------------------------------------------------------------------- #

# With the possibility of running as a batch job we do not want
# another build starting in this directory. Use a simple lock file
# to avoid multiple builds.
lock_file=.amanzi-build-lock
if [ -e $lock_file ]; then
    echo "Lock file ($lock_file) exists. Will not run until it is removed."
    exit_now "Fail to run with lock file present" 50
fi

# Function to remove the lock and set a trap to remove file.
# This trap will remove the lock file for any exit after this command.  
# Any exit called before this trap command will NOT remove the lock file.
function cleanup_lock()
{
    if [ -e $lock_file ] ; then
	echo "Removing the $lock_file"
	rm $lock_file
    fi
}
trap cleanup_lock 0 1 2 15

# Create the lock file 
echo "Creating lock file"
touch $lock_file

# ---------------------------------------------------------------------------- #
# Prep before the build                 
# ---------------------------------------------------------------------------- #

# Check structured and unstructured
if [ ${enable_structured} == "FALSE" -a ${enable_unstructured} == "FALSE" ]; then
    echo "Please specify structured or unstructured mesh capability"
    amanzi_build_usage
fi

# Submit as a batch job
if [ ${batch_build} == "TRUE" ]; then
    write_dflt_batch_file ${batch_build_file}
    submit_batch_script ${batch_build_file}
    echo "Do NOT delete or alter ${batch_build_file}"
    exit_now "Montior this directory for PBS log files" 0
fi

# Load modules
. $ASCEM_HOME/modules/init/ascem.bashrc ${prg_env}

# Define the init file 
if [ -z $init_file ]; then
    default_init_file
fi

# Check init file existence
if [ ! -e $init_file ]; then
    exit_now "CMake initialize file ${init_file} does not exist" 20
fi

# Update the repository branch
update_repository_branch ${amanzi_source_dir} ${amanzi_branch}

# Define the production install directory
if [[ "${amanzi_install_dir}" == "production" ]]; then
    echo "Running a production install"
    amanzi_install_dir="`production_install_dir`"
    amanzi_production=TRUE
fi

# Clobber old files
if [[ "$clobber" == "TRUE" ]]; then
    if [[ -e Makefile ]]; then
        echo "Cleaning build directory"
	make clean
    fi
    if [[ -e CMakeCache.txt ]]; then
	echo "Removing previous CMakeCache.txt file"
	rm CMakeCache.txt
    fi
fi

# Configure                
echo "Configure Amanzi in ${amanzi_source_dir}"
cmake \
    -C ${init_file} \
    -D CMAKE_INSTALL_PREFIX:FILEPATH=${amanzi_install_dir} \
    -D BUILD_SHARED_LIBS:BOOL=FALSE \
    -D TESTS_REQUIRE_MPIEXEC:BOOL=TRUE \
    -D TESTS_REQUIRE_FULLPATH:BOOL=TRUE \
    -D PREFER_STATIC_LIBRARIES:BOOL=TRUE \
    -D BUILD_STATIC_EXECUTABLES:BOOL=TRUE \
    -D ENABLE_Structured:BOOL=${enable_structured} \
    -D ENABLE_Unstructured:BOOL=${enable_unstructured} \
    -D ENABLE_HYPRE:BOOL=${amanzi_enable_hypre} \
    ${amanzi_source_dir}

if [ $? -ne 0 ]; then
    exit_now "Failed to configure Amanzi" 20
fi


# ---------------------------------------------------------------------------- #
# Build                 
# ---------------------------------------------------------------------------- #
make -j ${make_np}

if [ $? -ne 0 ]; then
    exit_now "Failed to build Amanzi" 20
fi

# ---------------------------------------------------------------------------- #
# Clean-up after the build                 
# ---------------------------------------------------------------------------- #

# Run test suite
if [[ "${amanzi_production}" == "TRUE" ]]; then
    echo "Now running test suite"
    id="`mercurial_global_id ${amanzi_source_dir}`"
    branch="`mercurial_branch ${amanzi_source_dir}`"
    test_results="${amanzi_install_dir}/bin/test-results.${id}.${branch}"
    echo "Results will be saved in ${test_results}"
    ctest --output-on-failure --output-log ${test_results}

    if [[ $? -ne 0 ]]; then
	echo 'Build did not pass the test suite. Will not install.'
	exit_now
    fi
fi

# Install
make install

if [ $? -ne 0 ]; then
    exit_now "Failed to install Amanzi" 20
fi

# Rename the binary if this is a production install
if [[ "${amanzi_production}" == "TRUE" ]]; then

    save_dir=`pwd`

    cd ${amanzi_install_dir}/bin

    # Move the binary and create a link to amanzi
    old_binary=amanzi
    new_binary=`production_binary_name ${amanzi_source_dir}`
    echo "Moving ${old_binary} to ${new_binary}"
    mv ${old_binary} ${new_binary}
    echo "Creating soft link"
    ln -s ${new_binary} ${old_binary}

    cd ${save_dir}

    # Change the permissions
    echo "Changing permissions in ${amanzi_install_dir}"
    chgrp -R ${ASCEM_GROUP} ${amanzi_install_dir}
    chmod -R g+rwX ${amanzi_install_dir}
    chmod -R g+s ${amanzi_install_dir}

fi


exit_now "$0 successfully finished" 0

