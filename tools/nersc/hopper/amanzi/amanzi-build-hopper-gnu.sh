#!/bin/bash

# ########################################################################### #
# Amanzi Build Script                                                         #
# GNU compiler programming environment                                        #
# ########################################################################### #


# --------------------------------------------------------------------------- #
# Module File Setup
# --------------------------------------------------------------------------- #

# Initialize the module command
if [[ $MODULESHOME ]]; then
    . $MODULESHOME/init/bash
fi

# Unload ENV module...want a clean ENV area
compiler_loaded=false
if [[ $PE_ENV ]]; then
    compiler_loaded=`echo $PE_ENV | tr "[A-Z]" "[a-z]"`
fi

if [[ $compiler_loaded ]]; then
    module unload PrgEnv-${compiler_loaded}
fi

# Add ASCEM module files
module use /project/projectdirs/m1012/modulefiles/${NERSC_HOST}
module load AmanziEnv-gnu

# --------------------------------------------------------------------------- #
# Process the command line
# --------------------------------------------------------------------------- #

# Define the defaults

# Actions
AMANZI_HELP=0
AMANZI_CLOBBER=0
AMANZI_CONFIG=0
AMANZI_MAKE=0
AMANZI_TEST=0

# Paths
AMANZI_SRC_DIR=
AMANZI_INSTALL_DIR=$HOME/amanzi-install/${NERSC_HOST}
AMANZI_BUILD_DIR=

# Parameters
AMANZI_NP_MAKE=2

# CCSE configuration

# 1, 2, 3
AMANZI_SPACEDIM=2

# FLOAT, DOUBLE
AMANZI_PRECISION=DOUBLE

# AMANZI, COREREACT (chooses chemistry evolution module)
AMANZI_CHEMEVOL_PKG=AMANZI


print_help()
{
echo '
Usage: '"$0"' [options]
Options: [defaults in brackets after descriptions]
Configuration:
  -h                      print this message
  -b                      Clobber Amanzi build directory tree
  -c                      Configure Amanzi (i.e. run cmake)
  -m                      Build and install Amanzi
  -t                      Run Amanzi test suite
  -a                      Configure, build and run test suite

Parameters:
  -n NP                   Number of parallel make threads

Directory and file names:
  -d SCR_DIR              Amanzi source directory
                          [${AMANZI_SRC_DIR}]
  -i INSTALL_DIR          Installation directory INSTALL_DIR
                          [${AMANZI_INSTALL_DIR}]

Run script this script in the directory you want to build Amanzi.
Currently only the unstructured Amanzi code is enabled with this script.
'
  exit 10
}



while getopts "habcd:i:mn:t" flag
do
  case $flag in
    h) AMANZI_HELP=1;;  
    a) AMANZI_CONFIG=1; AMANZI_MAKE=1; AMANZI_TEST=1;;
    b) AMANZI_CLOBBER=1;;
    c) AMANZI_CONFIG=1;;
    m) AMANZI_MAKE=1;;
    n) AMANZI_NP_MAKE=${OPTARG};;
    t) AMANZI_TEST=1;;
    d) AMANZI_SRC_DIR=${OPTARG};; 
    i) AMANZI_INSTALL_DIR=${OPTARG};; 
  esac
done

# Print help if -h was used 
if [[ $AMANZI_HELP -eq 1 ]]; then
    print_help
fi

# Check options
if [[ -z "${AMANZI_SRC_DIR}" ]]; then
    echo "Please define an Amanzi source directory with teh -d option"
    print_help
fi

if [[ ! -d $AMANZI_SRC_DIR ]]; then
    echo "Amanzi source directory \"${AMANZI_SRC_DIR}\" does not exist"
    print_help
fi
AMANZI_BUILD_DIR=${AMANZI_SRC_DIR}/build
mkdir -p ${AMANZI_BUILD_DIR}

if [ $AMANZI_NP_MAKE = "$(echo $AMANZI_NP_MAKE} | awk '{print strtonum($0)}')" ]; then
    echo "Running parallel make with $AMANZI_NP_MAKE threads"
else
    echo "Invalid -n option must be a number"
    print_help
fi


# Print the status after command line parsed
echo "----------------------------------------------"
echo "Amanzi source directory ${AMANZI_SRC_DIR}"
echo "Amanzi install directory ${AMANZI_INSTALL_DIR}"
echo ""
echo "Configuration Options"
echo "AMNZI_CLOBBER=${AMANZI_CLOBBER}"
echo "AMNZI_CONFIG=${AMANZI_CONFIG}"
echo "AMNZI_MAKE=${AMANZI_MAKE}"
echo "AMNZI_TEST=${AMANZI_TEST}"
echo ""
echo "Parameters"
echo "Number of make threads=${AMANZI_NP_MAKE}"
echo "----------------------------------------------"



# --------------------------------------------------------------------------- #
# Execute actions
# --------------------------------------------------------------------------- #

# Clobber
if [[ $AMANZI_CLOBBER -eq 1 ]]; then
    rm -rf ${AMANZI_BUILD_DIR}
    cd ${AMANZI_SRC_DIR}
    if [[ -e Makefile ]]; then
        make clean
    fi
    # explicitly specify src/*.cmake so we do not accidently remove *.cmake files from tools!
    cd src
    rm -rf \
        CMakeCache.txt \
        *.cmake */*.cmake */*/*.cmake \
        CMakeFiles */CMakeFiles */*/CMakeFiles \
        Makefile */Makefile */*/Makefile \
        Testing \
        *~ */*~ */*/*~ \
        *.rej */*.rej */*/*.rej
fi

# Configure
if [[ $AMANZI_CONFIG -eq 1 ]]; then

    export CXX=`which CC`
    export CC=`which cc`
    export BOOST_ROOT=${BOOST_DIR}

    echo "Running cmake to configure ${AMANZI_SRC_DIR}"

    cd ${AMANZI_BUILD_DIR}

    cmake \
        -D ENABLE_Config_Report:BOOL=ON \
	-D CMAKE_INSTALL_PREFIX:FILEPATH=${AMANZI_INSTALL_DIR} \
	-D BUILD_SHARED_LIBS:BOOL=OFF \
	-D PREFER_STATIC_LIBRARIES:BOOL=ON \
	-D BUILD_STATIC_EXECUTABLES:BOOL=ON \
	-D CMAKE_C_COMPILER=`which cc` \
	-D CMAKE_CXX_COMPILER=`which CC` \
	-D CMAKE_Fortran_COMPILER=`which ftn` \
        -D MPI_DIR:FILEPATH=${MPICH_DIR} \
        -D MPI_EXEC:FILEPATH=aprun \
        -D MPI_EXEC_NUMPROCS_FLAG:STRING=-n \
        -D MPI_EXEC_MAX_NUMPROCS_FLAG:STRING=4 \
        -D ENABLE_TESTS:BOOL=ON \
	-D BOOST_DIR:FILEPATH=${BOOST_DIR} \
        -D UnitTest_DIR:FILEPATH=${UnitTest_DIR} \
        -D HDF5_DIR:FILEPATH=${HDF5_DIR} \
	-D ENABLE_ASCEMIO:BOOL=ON \
	-D ASCEMIO_DIR:FILEPATH=${ASCEMIO_DIR} \
        -D NetCDF_DIR:FILEPATH=${NETCDF_DIR} \
        -D ExodusII_DIR:FILEPATH=${EXODUSII_DIR} \
        -D ENABLE_MESH_MOAB:BOOL=OFF \
        -D ENABLE_MESH_MSTK:BOOL=OFF \
        -D MSTK_DIR:FILEPATH=${MSTK_DIR} \
	-D MSTK_LIBRARY_DIR:FILEPATH=${MSTK_DIR}/lib/x86_64_Linux \
        -D METIS_DIR:FILEPATH=${METIS_DIR} \
        -D ENABLE_STK_Mesh:BOOL=ON \
        -D Trilinos_DIR:FILEPATH=${TRILINOS_DIR} \
	-D ENABLE_Unstructured:Bool=ON \
	-D ENABLE_Structured:Bool=ON \
        -D CCSE_DIR:FILEPATH=${CCSE_DIR} \
	-D ENABLE_MPI:Bool=ON \
	-D ENABLE_OpenMP:Bool=ON \
        -D AMANZI_SPACEDIM:INT=${AMANZI_SPACEDIM} \
        -D AMANZI_CHEMEVOL_PKG:STRING=${AMANZI_CHEMEVOL_PKG} \
        -D AMANZI_PRECISION:STRING=${AMANZI_PRECISION} \
        ..
    
    if [ $? -ne 0 ]; then
	echo "Failed to configure"
        exit 1
    fi
fi

# Build
if [[ $AMANZI_MAKE -eq 1 ]]; then
    cd ${AMANZI_BUILD_DIR}
    make -j ${AMANZI_NP_MAKE}
    if [ $? -ne 0 ]; then
	echo "Failed to build"
        exit 1
    else
	make install
	if [ $? -ne 0 ]; then
	    echo "Failed to install Amanzi"
	    exit 1
	fi    
    fi
fi

# Test
if [[ $AMANZI_TEST -eq 1 ]]; then
    echo "WARNING some tests will fail"
    cd ${AMANZI_BUILD_DIR}
    ctest --timeout 60 --output-on-failure
    if [ $? -ne 0 ]; then
        exit 1
    fi
fi

