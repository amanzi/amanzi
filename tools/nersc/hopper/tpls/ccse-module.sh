#!/bin/bash

# configure CCSE
CCSE_VERSION=0.1.7
CCSE_VERBOSE=0
SPACEDIM=2
ENABLE_MPI=1
ENABLE_OpenMP=1
PARALLEL_NP=1

# need this stuff for the module file
# should pull this automatically
COMPILER_BRAND=gcc
COMPILER_NAME=gcc
COMPILER_VERSION=4.6.1

CCSE_PREFIX=/project/projectdirs/m1012/tpls/ccse
CCSE_INSTALL_DIR=${CCSE_PREFIX}/${CCSE_VERSION}/${NERSC_HOST}-${COMPILER_NAME}-${COMPILER_VERSION}

BASE_DIR=`pwd`

################################################################################
#
# ccse
#
################################################################################
function install_ccse {
    CCSE_DIR=${CCSE_PREFIX}/ccse-${CCSE_VERSION}
    CCSE_CONFIG=1
    CCSE_MAKE=1
    CCSE_TEST=1
    CCSE_INSTALL=1

    # 1, 2, 3
    CCSE_SPACEDIM=${SPACEDIM}

    # FLOAT, DOUBLE
    CCSE_PRECISION=DOUBLE

    CCSE_BUILD_DIR=${CCSE_DIR}/build-${CCSE_SPACEDIM}D
    rm -rf ${CCSE_BUILD_DIR}
    mkdir -p ${CCSE_BUILD_DIR}
    cd ${CCSE_BUILD_DIR}


    if [ ${CCSE_VERBOSE} -eq 1 ]; then
	VFLAG="--debug-output"
    fi

    export CC=`which cc`
    export CXX=`which CC`
    export FC=`which ftn`

    cmake \
	-D CMAKE_INSTALL_PREFIX:FILEPATH=${CCSE_INSTALL_DIR} \
	-D ENABLE_Config_Report:BOOL=ON \
        -D CMAKE_C_COMPILER:FILEPATH=$CC \
        -D CMAKE_CXX_COMPILER:FILEPATH=$CXX \
        -D CMAKE_Fortran_COMPILER:FILEPATH=$FC \
	-D ENABLE_MPI:BOOL=${ENABLE_MPI} \
	-D ENABLE_OpenMP:BOOL=${ENABLE_OpenMP} \
	-D ENABLE_TESTS:BOOL=ON \
	-D BL_SPACEDIM:INT=${CCSE_SPACEDIM} \
	-D BL_PRECISION:STRING="${CCSE_PRECISION}" \
	${VFLAG} \
	..

    if [ $? -ne 0 ]; then
	exit 1
    fi

    if [ ${CCSE_VERBOSE} -eq 1 ]; then
	VFLAG="VERBOSE=ON"
    fi
    make -j ${PARALLEL_NP} ${VFLAG}
    if [ $? -ne 0 ]; then
	exit 1
    fi

    #ctest --timeout 60 --output-on-failure
    if [ $? -ne 0 ]; then
	exit 1
    fi

    make install
    if [ $? -ne 0 ]; then
	exit 1
    fi
    
}

################################################################################
#
# write the module file
#
# copy this file to ${project_home}/modulefiles/${HOST}/ccse/${CCSE_VERSION}
#
################################################################################
function write_module_file {
    echo "Writing module file...."

    cat > modulefile-${CCSE_VERSION} <<EOF
#%Module###<-magic cookie ####################################################
##
## Amanzi CCSE ${CCSE_VERSION} Module
##

# Define conflicts here
conflict ccse

# Local Tcl variables
set module_name                  ccse
set version                      ${CCSE_VERSION}
set ascem_home                   /project/projectdirs/m1012

# Set the machine name
if { [ info exists env(NERSC_HOST) ] } {
    set machine \$env(NERSC_HOST)
} else {
    set machine not_found
}


# List of supported compilers
# Unfortunately if the default is loaded 
# the PE_ENV is not set when this runs and 
# therefore we must define the version here.
set valid_compiler_list      [list gnu]
set default_compiler         gnu
set default_compiler_version ${COMPILER_VERSION}

#ModulesHelp -- output from module help module_name/version
proc ModulesHelp { } {

    global module_name
    global version
    global ascem_home
    global valid_compiler_list

    puts stderr ""
    puts stderr "==================================================================="
    puts stderr "\$module_name \$version"
    puts stderr "Location: \$ascem_home/tpls/\$module_name/\$version"
    puts stderr "==================================================================="
    puts stderr "\$module_name \$version built for the Amanzi code project"
    puts stderr "Built with parallel support"
    puts stderr "Supported Compiler Environments: \$valid_compiler_list"
    puts stderr ""
    puts stderr "===================================================================\n"

}

module-whatis   "LBNL CCSE BoxLib library."



# Determine the compiler, use the default if a PrgEnv-* is not loaded
set compiler          not_found
if { [ info exists env(PE_ENV) ] } {
    set compiler_lc [ string tolower \$env(PE_ENV) ]
    if { [ lsearch -exact \$valid_compiler_list \$compiler_lc ] >= 0 } {
	set compiler \$compiler_lc
    } else {
	puts stderr "\${env(PE_ENV)} is not a supported compiler environment"
	puts stderr "Please load a supported compiler environment"
	puts stderr "Supoorted compiler environments: \$valid_compiler_list"
	break 
	exit 1
    }
} else {
    set compiler \$default_compiler
    module load PrgEnv-\$compiler
}

if { [ string equal \$compiler "not_found" ] } {
    puts stderr "Failed to determine compiler. Can not define install prefix."
    break 
    exit 1
}

# Determine the compiler version
set compiler_version         0
if { [ info exists env(PE_ENV) ] } {
    set compiler_version_env_name  \${env(PE_ENV)}_VERSION
    if { [ info exists env(\$compiler_version_env_name) ] } {
	set compiler_version [ expr { \$env(\$compiler_version_env_name) } ]
    } else {
	puts stderr "Although PE_ENV=\$env(PE_ENV) is set, can not determine compiler version"
    }
} else {
    if { [ string equal \$compiler \$default_compiler ] } {
	set compiler_version \$default_compiler_version
    } else {
	puts stderr "PE_ENV is not set and \$compiler is not the default compiler"
    }
}

if { "\$compiler_version" == "0" } {
    puts stderr "Failed to determine the compiler version, can not define install prefix"
    break 
    exit 1
}
    

# Determine the build type
# Build Type machine-compiler-complier_version
# Now set the build type based on machine, compiler and compiler version
# GNU is a special case. We use gcc instead of gnu in the build type definition
if { [ string equal "\$compiler" "gnu" ] } {
    set build_type \$machine-gcc-\$compiler_version
} else {
    set build_type \$machine-\$compiler-\$compiler_version
}

# Now set the root install path
set install_prefix              \$ascem_home/tpls/\$module_name/\$version/\$build_type

# Complain if the prefix directory does not exist
if { [ file exists \$install_prefix ] } {
} else {
    puts stderr "Install prefix \$install_prefix does not exist or user has insufficient privileges to search"
    break
    exit 1
}

# Define the lib, include, man paths based on the install prefix
set lib_path                        \$install_prefix/lib
set include_path                    \$install_prefix/include

# Environment variables
setenv CCSE_VERSION             \$version
setenv CCSE_DIR                 \$install_prefix

# *PATH environment variables
prepend-path LD_LIBRARY_PATH        \$lib_path

EOF
}

################################################################################
#
# main
#
################################################################################

# module setup
. ${MODULESHOME}/init/bash

AMANZI_HOME=/project/projectdirs/m1012
AMANZI_TPL_DIR=${AMANZI_HOME}/tpls

module use ${AMANZI_HOME}/modulefiles/${NERSC_HOST}
module load PrgEnv-gnu
module load cmake

# write setup info to file:
DUMP_FILE=log-ccse-setup-`date +"%F_%k-%M"`
module list &> $DUMP_FILE
echo  >> $DUMP_FILE
echo "CCSE_VERSION = ${CCSE_VERSION}" >> $DUMP_FILE
echo "CCSE_PREFIX = ${CCSE_PREFIX}" >> $DUMP_FILE
echo "CCSE_INSTALL_DIR = ${CCSE_INSTALL_DIR}" >> $DUMP_FILE

install_ccse
cd ${BASE_DIR}
write_module_file