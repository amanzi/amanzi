#!/bin/sh
################################################################################
#                                                                              #
# Amanzi TPL Build Script                                                      #
#                                                                              #
################################################################################

# ---------------------------------------------------------------------------- #
# Initialize                  
# ---------------------------------------------------------------------------- #

# Variables to improve code readability
TRUE=1
FALSE=0

# Initialize the module command
if [[ $MODULESHOME ]]; then
    . $MODULESHOME/init/bash
fi

# Amanzi source directory
amanzi_source_dir=`pwd`

# Program environment
prg_env=gnu

# TPL install prefix
tpl_install_prefix=

# Build type
tpl_build_type=Debug

# Number of make jobs
make_np=1

# ---------------------------------------------------------------------------- #
# Environment Variables                  
# ---------------------------------------------------------------------------- #
export ASCEM_HOME=/project/projectdirs/m1012
export ASCEM_TPL_INSTALL_DIR=${ASCEM_HOME}/tpls/install/${NERSC_HOST}

# ---------------------------------------------------------------------------- #
# Useful functions                  
# ---------------------------------------------------------------------------- #
function exit_now()
{
    echo $1
    exit $2
}
function print_usage()
{
echo '
Usage: '"$0"' [options] DIR
where Amanzi source is located in DIR ['"${amanzi_source_dir}"']
Options: [Defaults in brackets]

    --help                  print this message
    --prg-env=STRING        load PrgEnv module STRING ['"${prg_env}"']        
    --parallel=n            number of parallel make threads ['"${make_np}"']
    --install=DIR           install Amanzi TPL in DIR
    --debug                 build debug libraries
    --opt                   build optimized libraries
'
  exit 10
}
function pe_env_tolower()
{
    echo ${PE_ENV} | tr "[A-Z]" "[a-z]"
}
function get_compiler_version()
{
    version_var=${PE_ENV}_VERSION
    echo ${!version_var}
}
get_mpi_version()
{
    echo ${CRAY_MPICH2_VER}
}
load_prg_env()
{
    compiler_loaded=false
    if [[ $PE_ENV ]]; then
	compiler_loaded=`echo $PE_ENV | tr "[A-Z]" "[a-z]"`
    fi

    if [[ "$compiler_loaded" != "false" ]]; then
	if [[ "${compiler_loaded}" != $1 ]] ; then
	  echo "Swapping ${compiler_loaded} with $1"
	  module swap PrgEnv-${compiler_loaded} PrgEnv-$1
        fi  
    else
	module load PrgEnv-$1
    fi

    echo "Loading PrgEnv-$1"
    module load PrgEnv-$1
    if [[ $1 == "gnu" ]]; then
	module swap gcc gcc/4.6.1
    fi
}
default_install_prefix()
{
   mpi_version=`get_mpi_version` 
   comp_version=`get_compiler_version`

   echo "${ASCEM_TPL_INSTALL_DIR}/mpich2-${mpi_version}-${prg_env}-${comp_version}/${tpl_build_type}"

}

function parse_option_with_equal()
{
  a=$1
  opt_name=$2

  if echo $a | grep "^--${opt_name}=" > /dev/null 2> /dev/null; then
    echo $a | sed "s/^--${opt_name}=//"
  fi

}

function parse_argv()
{
  argv=( "$@" )
  last=$(( ${#argv[@]} - 1 ))
  i=0
  while [ $i -le ${last} ]
  do
    opt=${argv[$i]}
    #echo "i: ${i} opt=$opt last: $last"

    # Skip the case search is this is a directory
    if [[ -d $opt ]]; then
	amanzi_source_dir=$opt
        i=$[$i+1]
	continue
    fi

    case ${opt} in

      -h|--h|--help)
                print_usage
		exit_now 0
		;;

      --prg-env=*)
                 prg_env=`parse_option_with_equal ${opt} 'prg-env'`
		 ;;

      --install=*)
                 tpl_install_prefix=`parse_option_with_equal ${opt} 'install'`
		 ;;

      --parallel=[0-9]*)
                 parallel_jobs=`parse_option_with_equal ${opt} 'parallel'`
		 ;;

      --debug)
                 tpl_build_type=Debug
		 ;;

      --opt)
                 tpl_build_type=Release
		 ;;

          *)
		   exit "'${opt}' is an unknown option" 20
		   ;;
      esac

      i=$[$i+1]
  done

}

#
# End functions
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
#  Main 
# ---------------------------------------------------------------------------- #

# Parse the command line arguments
array=( "$@" )
parse_argv "${array[@]}"

# Load the programming module
#load_prg_env $prg_env
# Use the ASCEM init file to load the correct modules
. $ASCEM_HOME/modules/init/ascem.bashrc $prg_env

# Check the install prefix
if [[ -z "${tpl_install_prefix}" ]]; then
    tpl_install_prefix=`default_install_prefix`
fi

# TPL project root
tpl_project_root=${amanzi_source_dir}/config/SuperBuild

# Define the libsci config file
libsci_file=${tpl_project_root}/include/trilinos-blas-libsci-${prg_env}.cmake

echo ""
echo "--------------------------------------------------------------------------------"
echo " Building Amanzi TPL"
echo " Build Type:             ${tpl_build_type}"
echo " TPL CMake Project Root: ${tpl_project_root}"
echo " Trilinos Config File:   ${libsci_file}"
echo " Install Location:       ${tpl_install_prefix}"
echo ""
echo -n " C Compiler:             "
which cc
echo -n " C++ Compiler:           "
which CC
echo -n " Fortran Compiler:       "
which ftn
echo "--------------------------------------------------------------------------------"
echo ""

# Module list
echo ""
echo "--------------------------------------------------------------------------------"
module list
echo "--------------------------------------------------------------------------------"
echo ""

# Configure (CMake) 
cmake \
    -D CMAKE_C_COMPILER:FILEPATH=`which cc` \
    -D CMAKE_CXX_COMPILER:FILEPATH=`which CC` \
    -D CMAKE_Fortran_COMPILER:FILEPATH=`which ftn` \
    -D TPL_INSTALL_PREFIX:PATH=${tpl_install_prefix} \
    -D CMAKE_BUILD_TYPE:STRING=${tpl_build_type} \
    -D MPI_EXEC:STRING='/usr/bin/aprun' \
    -D ENABLE_HYPRE:BOOL=TRUE \
    -D Trilinos_Build_Config_File:FILEPATH=${libsci_file} \
    ${tpl_project_root}

if [ $? -ne 0 ] ; then
    echo "CMake configure failed"
    exit $?
fi

# Build
make -j ${make_np}

if [ $? -ne 0 ] ; then
    echo "Make step failed"
    exit $?
fi

# Install
make install


if [ $? -ne 0 ] ; then
    echo "Install failed"
fi

exit $?
