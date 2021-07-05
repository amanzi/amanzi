#!/bin/bash

# ############################################################################ #
#
#  Amanzi Third Party Configurations
#
# ############################################################################ #
#
# This script is intended to be a temporary document for third party library
# configurations. It is a placeholder for Python scripts that are
# currently under development. 
#
#
# Notes:
# 
# () User should copy this script to a private copy and modify the copy. It is
#    not intended to be a one-script-for-all-systems tool. It is a placeholder
#    for Python based tools that are currently under development.
#
# () The configuration options defined in each *_config function are a common
#    set of choices that have been successful on multiple machines, but it does
#    guarantee that these particular choices will work on any given system.
#    Please do not push changes to this script that contain personal
#    configuration options. 
#   
# () This script expects the distribution files for each of the TPLs in the same
#     directory. This path is defined by the DOWNLOAD_PREFIX variable.
#
# () All of the TPLs will be built and installed in the same directory. The 
#    variable PREFIX defines the directory for the builds and install 
#    directories.
#
# () To skip a TPL, set the *_PREFIX variable to value other then PREFIX
#
# () Amanzi does not support non-MPI builds. The environment variables CC and
#    CXX are set to the MPI compiler wrappers mpicc and mpicxx. The user is
#    required to provide the location of the MPI installation with the 
#    MPI_PREFIX variable.
#
# () The variable PLATFORM gates configurations for particular OS types.
#    The recognized values are: MacOS, RedHat, Fedora, Ubuntu and Fedora

# System Libraries

# CMake (Need version 2.8.2 or higher)
CMAKE=cmake

# Number of parallel make jobs (make -j ${MAKE_NP})
MAKE_NP=4

#
# MPI (Amanzi does not support non-MPI builds)
# Set
MPI_PREFIX=${MPI_ROOT}
MPICC=${MPI_PREFIX}/bin/mpicc
MPICXX=${MPI_PREFIX}/bin/mpicxx
MPIEXEC=${MPI_PREFIX}/bin/mpiexec
MPIEXEC_NUMPROCS_FLAG=-np
MPIEXEC_MAX_NUMPROCS=4

# Platform (Platforms:MacOS,RedHat,Ubuntu,Fedora) us this to gate 
# Use this variable to gate configuration requirements for particular
# platforms.
PLATFORM=

# 
# Define the prefix directories
#
if [ !  ${PREFIX} ]; then
	PREFIX=`pwd`
fi
SOURCE_PREFIX=${PREFIX}/src
BUILD_PREFIX=${PREFIX}/build

#
# Define the download directory
# Location of tar and zip source files
DOWNLOAD_DIRECTORY=$PWD

#
# Amanzi Build Script
#
# This script will generate a another bash shell
# that will build amanzi using the TPL
# locations defined in this script
#
BUILD_AMANZI_SCRIPT=0
AMANZI_BUILD_SCRIPT=amanzi-build.sh

#
# TPL
#
# Each TPL will have 3 variables defined
# BUILD_XXX, XXX_PREFIX and XXX_VERSION
# The build variable indicates if the package should
# be built by this script. The PREFIX variable
# defines the installation prefix. If the script
# does not build the TPL, the user MUST provide the
# installation prefix to build dependent packages.
# The version variable controls the name of the
# tarfile or zip source file. 

#
# Boost
#
BUILD_BOOST=0
BOOST_PREFIX=${PREFIX}
BOOST_VERSION=1.48.0

#
# UnitTest
#
BUILD_UNITTEST=0
UNITTEST_PREFIX=${PREFIX}
UNITTEST_VERSION=1.4

#
# ZLIB 
#
BUILD_ZLIB=0
ZLIB_PREFIX=${PREFIX}
ZLIB_VERSION=1.2.6

#
# cURL
#
BUILD_CURL=0
CURL_PREFIX=${PREFIX}
CURL_VERSION=7.24.0

#
# HDF5
#
BUILD_HDF5=0
HDF5_PREFIX=${PREFIX}
HDF5_VERSION=1.8.8

#
# ASCEM-IO
#
BUILD_ASCEMIO=0
ASCEMIO_PREFIX=${PREFIX}
ASCEMIO_VERSION=1.2

#
# netCDF
#
BUILD_NETCDF=0
NETCDF_PREFIX=${PREFIX}
NETCDF_VERSION=4.1.1

#
# ExodusII
#
BUILD_EXODUS=0
EXODUS_PREFIX=${PREFIX}
EXODUS_VERSION=5.14

#
# METIS
#
BUILD_METIS=0
METIS_PREFIX=${PREFIX}
METIS_VERSION=4.0.3

#
# MSTK
#
BUILD_MSTK=0
MSTK_PREFIX=${PREFIX}
MSTK_VERSION=1.83

#
# MOAB
# MOAB_VERSION=moab_svn will download the latest version
BUILD_MOAB=0
MOAB_PREFIX=${PREFIX}
MOAB_VERSION=r4276

#
# Trilinos
#
BUILD_TRILINOS=0
TRILINOS_PREFIX=${PREFIX}
TRILINOS_VERSION=10.6.4

################################################################################
#
# Very little below this point should need to be changed....
#
#################################################################################

SCRIPT_DIR=$PWD

TAR_FLAGS=zxf

################################################################################
#
# FUNCTIONS
#
#################################################################################

# ################################### #
# help_message                        #
# ################################### #
function help_message {
    echo "Usage: ${0##*/} [OPTION]
    Default action is to not build any TPL and generate the Amanzi build script
    -h              Print this help message and exit
    -d <path>       Location of tar and zip source files (DOWNLOAD_DIRECTORY)            

    TPL build flags
    -a              Build all TPLs
    -A              Build ASCEM-IO
    -b              Build Boost
    -e              Build Exodus
    -H              Build HDF5
    -k              Build MSTK
    -m              Build MOAB
    -n              Build netCDF
    -s              Build METIS
    -t              Build Trilinos
    -u              Build UnitTest
    -w              Build Curl
    -z              Build ZLIB

    User is required to provide location of the MPI installation by defining the
    variable MPI_PREFIX in this script.

    TPL built by this script must locate all tar and zip source files in one 
    directory. This directory is defined by the variable 
    DOWNLOAD_DIRECTORY=${DOWNLOAD_DIRECTORY}.
    "

    exit 0

}

# ################################### #
# variable_is_set                     #
# ################################### #
function variable_is_set() {
echo "In function variable_is_set $1"
    [[ ! ${!1} && ${!1-_} ]] 
}
# ################################### #
# check_platform                      #
# ################################### #
function check_platform {

    if [ $PLATFORM ]; then
	echo "User defined platform ${PLATFORM}"
    else
	OS=`uname -s`
	if [ $OS == Darwin ]; then
	    PLATFORM=MacOS
	elif [ $OS == Linux ]; then
	    if [ -e /etc/lsb-release ]; then
	        # probably ubuntu.... may need to grep the file?
	        PLATFORM=Ubuntu
	    elif [ -e /etc/fedora-release ]; then
	        PLATFORM=Fedora
	    elif [ -e /etc/redhat-release ]; then
	        PLATFORM=RedHat
            else
		echo "Unknown Linux system"
		PLATFORM=Linux
            fi
	fi
	echo "Running on ${PLATFORM} platform"
    fi
}

# ################################### #
# check_mpi                           #
# ################################### #
function check_mpi {
  
    if [[  $MPI_PREFIX  &&  -e $MPI_PREFIX ]]; then
	echo "MPI installation located at ${MPI_PREFIX}"
	if [[ $MPICC && -e $MPICC ]]; then
	    echo "C Compiler info"
	    ${MPICC} -show
	else
	    echo "Can not locate mpicc"
	    help_message
	fi
	if [[ $MPICXX  &&   -e ${MPICXX} ]]; then
	    echo "C++ Compiler info"
	    ${MPICXX} -show
	else
	    echo "Can not locate mpicxx C++ compiler"
	    help_message
	fi
    else
	echo "MPI installation is not defined"
	help_message
    fi

}

# ################################### #
# check_download                      #
# ################################### #
function check_download {

    if [[ $DOWNLOAD_DIRECTORY && -e $DOWNLOAD_DIRECTORY ]]; then
	echo "Search source files in $DOWNLOAD_DIRECTORY"
    else
	echo "Must define a download directory"
	help_message
    fi
}

# ################################### #
# check_prefix                        #
# ################################### #
function check_prefix {

    if [[ $PREFIX ]]; then

	if [[ ! -e $PREFIX ]]; then
	    mkdir -p $PREFIX
	fi
    else
	echo "Must define a TPL PREFIX"
	help_message
    fi

    if [[ $SOURCE_PREFIX ]]; then

	if [[ ! -e $SOURCE_PREFIX ]]; then
	    mkdir -p $SOURCE_PREFIX
	fi
    else
	echo "Must define a TPL SOURCE_PREFIX"
	help_message
    fi

    if [[ $BUILD_PREFIX ]]; then

	if [[ ! -e $BUILD_PREFIX ]]; then
	    mkdir -p $BUILD_PREFIX
	fi
    else
	echo "Must define a TPL BUILD_PREFIX"
	help_message
    fi

}

# ################################### #
# print_build_status                  #
# ################################### #
function print_build_status {

    echo "BUILD_BOOST=$BUILD_BOOST        BOOST_PREFIX=$BOOST_PREFIX"
    echo "BUILD_CURL=$BUILD_CURL          CURL_PREFIX=$CURL_PREFIX"
    echo "BUILD_ZLIB=$BUILD_ZLIB          ZLIB_PREFIX=$ZLIB_PREFIX"
    echo "BUILD_HDF5=$BUILD_HDF5          HDF5_PREFIX=$HDF5_PREFIX"
    echo "BUILD_ASCEMIO=$BUILD_ASCEMIO    ASCEMIO_PREFIX=$ASCEMIO_PREFIX"
    echo "BUILD_NETCDF=$BUILD_NETCDF      NETCDF_PREFIX=$NETCDF_PREFIX"
    echo "BUILD_EXODUS=$BUILD_EXODUS      EXODUS_PREFIX=$PREFIX"
    echo "BUILD_MOAB=$BUILD_MOAB          MOAB_PREFIX=$PREFIX"
    echo "BUILD_METIS=$BUILD_METIS        METIS_PREFIX=$PREFIX"
    echo "BUILD_MSTK=$BUILD_MSTK          MSTK_PREFIX=$PREFIX"
    echo "BUILD_UNITTEST=$BUILD_UNITTEST  UNITTEST_PREFIX=$UNITTEST_PREFIX"
    echo "BUILD_TRILINOS=$BUILD_TRILINOS  TRILINOS_PREFIX=$PREFIX"

}

################################################################################
#
# curl
#
################################################################################
function build_curl {
    if [ ${CURL_PREFIX} == ${PREFIX} ]; then
        CURL_DIR=${PREFIX}/curl/curl-${CURL_VERSION}
        rm -rf ${CURL_DIR}
        mkdir -p ${PREFIX}/curl
        tar ${TAR_FLAGS} ${DOWNLOAD_DIRECTORY}/curl-${CURL_VERSION}.tar.gz -C ${PREFIX}/curl
        cd ${CURL_DIR}
        ./configure --prefix=${PREFIX} --enable-static
        if [ $? -ne 0 ]; then
            exit 
        fi
        make -j${MAKE_NP} all
        if [ $? -ne 0 ]; then
            exit 
        fi
        make install
        if [ $? -ne 0 ]; then
            exit 
        fi
    else
        echo Using curl from ${CURL_PREFIX}
    fi
}
################################################################################
#
# zlib
#
# Note: MacOSX: if zlib is built with mpicc rather than the system
# gcc, it will build .so files instead of .dylib files.... This cause
# link errors in trilinos and amanzi.
#
################################################################################
function build_zlib {
    if [ ${ZLIB_PREFIX} == ${PREFIX} ]; then
        ZLIB_DIR=${PREFIX}/zlib/zlib-${ZLIB_VERSION}
        rm -rf ${ZLIB_DIR}
        mkdir -p ${PREFIX}/zlib
	tarfile=${DOWNLOAD_DIRECTORY}/zlib-${ZLIB_VERSION}.tar.gz
        tar ${TAR_FLAGS}  ${tarfile} -C ${PREFIX}/zlib
        cd ${ZLIB_DIR}
        ./configure \
            --prefix=${PREFIX}

        if [ $? -ne 0 ]; then
            exit 
        fi
        # zlib doesn't play nice with parallel builds...?
        make all
        if [ $? -ne 0 ]; then
            exit 
        fi
        make install
        if [ $? -ne 0 ]; then
            exit 
        fi
    else
        echo Using zlib from ${ZLIB_PREFIX}
    fi
}


################################################################################
#
# unittest++
#
################################################################################
function build_unittest {
    if [ ${UNITTEST_PREFIX} == ${PREFIX} ]; then
        UNITTEST_DIR=${PREFIX}/unittest/UnitTest++
        rm -rf ${UNITTEST_DIR}
        mkdir -p ${PREFIX}/unittest
        unzip ${DOWNLOAD_DIRECTORY}/unittest-cpp-${UNITTEST_VERSION}.zip -d ${PREFIX}/unittest

        cd ${UNITTEST_DIR}
        make -j${MAKE_NP} all
        if [ $? -ne 0 ]; then
            exit 
        fi
        # ugh... manual install
        echo Copying UnitTest++ files to ${UNITTEST_PREFIX}
	if [ ! -e ${UNITTEST_PREFIX}/lib ]; then mkdir -p ${UNITTEST_PREFIX}/lib; fi
        cp libUnitTest++.a ${UNITTEST_PREFIX}/lib
	if [ ! -e ${UNITTEST_PREFIX}/include ]; then mkdir -p ${UNITTEST_PREFIX}/include; fi
        mkdir -p ${UNITTEST_PREFIX}/include/unittest++/Posix
        cp src/*.h ${UNITTEST_PREFIX}/include/unittest++
        cp src/Posix/*.h ${UNITTEST_PREFIX}/include/unittest++/Posix
        if [ $? -ne 0 ]; then
            exit 
        fi
    else
        echo Using unittest++ from ${UNITTEST_PREFIX}
    fi
}
################################################################################
#
# Boost
#
################################################################################
function build_boost {
    if [ ${BOOST_PREFIX} == ${PREFIX} ]; then
        file_version=`echo $BOOST_VERSION | sed "s:\.:_:g"`
	BOOST_DIR=${PREFIX}/boost/boost_${file_version}
        rm -rf ${BOOST_DIR}
        mkdir -p ${PREFIX}/boost
	tar ${TAR_FLAGS} ${DOWNLOAD_DIRECTORY}/boost_${file_version}.tar.gz -C ${PREFIX}/boost

	cd $BOOST_DIR
	./bootstrap.sh --prefix=${PREFIX}
	echo Boost bootstrap complete
	./bjam --prefix=$PREFIX install

    else
	echo Using boost form ${BOOST_PREFIX} 
    fi

}
################################################################################
#
# hdf5
#
# parallel is not compatible with c++ 
# (http://www.hdfgroup.org/hdf5-quest.html#p5thread)
# must disable the C++ interfaces since XDMF and ASCEM IO
# need parallel HDF5
#
################################################################################
function build_hdf5 {
    if [ ${HDF5_PREFIX} == ${PREFIX} ]; then
        HDF5_DIR=${PREFIX}/hdf5/hdf5-${HDF5_VERSION}
        rm -rf ${HDF5_DIR}
        mkdir -p ${PREFIX}/hdf5
        tar ${TAR_FLAGS} ${DOWNLOAD_DIRECTORY}/hdf5-${HDF5_VERSION}.tar.gz -C ${PREFIX}/hdf5

        cd ${HDF5_DIR}

        ./configure --prefix=${PREFIX} \
            --disable-fortran \
	        --disable-cxx \
            --enable-production \
            --enable-largefile \
            --enable-parallel 

        if [ $? -ne 0 ]; then
            exit 
        fi
        make -j ${MAKE_NP}
        if [ $? -ne 0 ]; then
            exit 
        fi
        if [ $BUILD_CHECK -eq 1]; then
            make check 
            if [ $? -ne 0 ]; then
                exit
            fi
        fi
        make install
    else
        echo Using hdf5 from ${HDF5_PREFIX}
    fi   
}
################################################################################
#
# ASCEM-IO
# 
# http://ascem-io.secure-water.org/
#
# Parallel IO library Requires HDF5
#
################################################################################
function download_ascemio {

   ASCEMIO_SVN_ROOT="//ascem-io.secure-water.org/ascem-io"
   if [ ${ASCEMIO_VERSION} == "trunk" ]
   then
	   ASCEMIO_DOWNLOAD_LOCATION="${ASCEMIO_SVN_ROOT}/trunk"
   else
	   ASCEMIO_DOWNLOAD_LOCATION="${ASCEMIO_SVN_ROOT}/releases/${ASCEMIO_VERSION}"
   fi

   echo "Download ASCEM-IO from ${ASCEMIO_DOWNLOAD_LOCATION}"

   svn co http:${ASCEMIO_DOWNLOAD_LOCATION} --username amanzi_dev --password gr0undw@t3r

   if [ $? -ne 0 ]; then
	   echo "Failed to download ASCEMIO"
	   exit
   fi

}
    
function build_ascemio {
    echo "Build ASCEMIO"

    if [ ${ASCEMIO_PREFIX} == ${PREFIX} ]; then
		SAVE_DIR=`pwd`
		ASCEMIO_DIR=${PREFIX}/ascem-io
		if [ ! -e ${ASCEMIO_DIR} ]; then
		   	mkdir -p ${ASCEMIO_DIR}
		fi
		ASCEMIO_SRC_DIR=${ASCEMIO_DIR}/src
		if [ ! -e ${ASCEMIO_SRC_DIR} ]; then
		   	mkdir -p ${ASCEMIO_SRC_DIR}
		fi
		cd ${ASCEMIO_SRC_DIR}
        download_ascemio
		cd ${ASCEMIO_VERSION}/src
	
		export CC=${MPICC}
		export CXX=${MPICXX}
		echo "Building ASCEM-IO against HDF5 located in ${HDF5_PREFIX}"
        make HDF5_INCLUDE_DIR=${HDF5_PREFIX}/include
		if [ $? -ne 0 ]; then
			echo "Failed to build ASCEM-IO"
			exit
		fi
		make ASCEMIO_INSTALL_DIR=${ASCEMIO_PREFIX} install
		if [ $? -ne 0 ]; then
			echo "Failed tp install ASCEM-IO package"
			exit
		fi
		cd ${SAVE_DIR}
	fi
    
}
################################################################################
#
# netcdf
#
################################################################################
function build_netcdf {
    if [ ${NETCDF_PREFIX} == ${PREFIX} ]; then
        NETCDF_DIR=${PREFIX}/netcdf/netcdf-${NETCDF_VERSION}
        rm -rf ${NETCDF_DIR}
        mkdir -p ${PREFIX}/netcdf
        tar ${TAR_FLAGS} ${DOWNLOAD_DIRECTORY}/netcdf-${NETCDF_VERSION}.tar.gz -C ${PREFIX}/netcdf
        
        cd ${NETCDF_DIR}
       
	# Must modify netCDF to build ExodusII
	echo "Patching the netCDF source"
        perl -w -i -p -e 's@#define NC_MAX_DIMS[\s]+[\d]+@#define NC_MAX_DIMS 65536@' libsrc/netcdf.h libsrc4/netcdf.h libsrc4/netcdf_base.h
        perl -w -i -p -e 's@#define NC_MAX_VARS[\s]+[\d]+@#define NC_MAX_VARS 524288@' libsrc/netcdf.h libsrc4/netcdf.h libsrc4/netcdf_base.h
        perl -w -i -p -e 's@#define NC_MAX_VAR_DIMS[\s]+NC_MAX_DIMS@#define NC_MAX_VAR_DIMS 8@' libsrc/netcdf.h libsrc4/netcdf.h libsrc4/netcdf_base.h

#        grep -e NC_MAX_DIMS libsrc4/netcdf.h
#        grep -e NC_MAX_VARS libsrc4/netcdf.h
#        grep -e NC_MAX_VAR_DIMS libsrc4/netcdf.h

	USE_NETCDF4='--disable-netcdf-4'
	if [ $PLATFORM == Ubuntu ]; then
	    USE_NETCDF4='--enable-netcdf-4'
	fi

        ./configure --prefix=${PREFIX} \
            --disable-fortran \
            --disable-f90 \
            --disable-f77 \
            --disable-fortran-compiler-check \
            ${USE_NETCDF4} \
            --disable-cxx-4 \
            --disable-dap \
            --with-mpi=${MPI_PREFIX} \
            --with-hdf5=${HDF5_PREFIX} \
            --with-zlib=${ZLIB_PREFIX}
        if [ $? -ne 0 ]; then
            exit 
        fi
        make -j ${MAKE_NP}
        if [ $? -ne 0 ]; then
            exit 
        fi
        #make check 
        if [ $? -ne 0 ]; then
            exit
        fi
        make install
    else
        echo Using netcdf from ${NETCDF_PREFIX}
    fi
}
################################################################################
#
# exodus
#
################################################################################
function build_exodus_cmake {
    # download source from: http://sourceforge.net/projects/exodusii/
    export NETCDF_DIR=${NETCDF_PREFIX}

    EXODUS_DIR=${PREFIX}/exodusii/exodusii-${EXODUS_VERSION}
    rm -rf ${EXODUS_DIR}
    mkdir -p ${PREFIX}/exodusii
    tar ${TAR_FLAGS} ${DOWNLOAD_DIRECTORY}/exodusii-${EXODUS_VERSION}.tar.gz -C ${PREFIX}/exodusii
    mkdir -p ${EXODUS_DIR}/build


    cd ${EXODUS_DIR}
    # fedora: possible error in next line \{NETCDF_LIBRARY\} should be {NETCDF_LIBRARY} is this failing on osx and ubuntu?
    # Fails on red hat as well, not clear why we are substituting here.
    #perl -w -i -p -e "s@TARGET_LINK_LIBRARIES\(exoIIv2c \$\{NETCDF_LIBRARY\}/libnetcdf\.a\)@TARGET_LINK_LIBRARIES\(exoIIv2c \${NETCDF_LIBRARY}/libnetcdf\.a \${NETCDF_LIBRARY}/libhdf5_hl\.a \${NETCDF_LIBRARY}/libhdf5\.a\)@" cbind/CMakeLists.txt

    cd ${EXODUS_DIR}/build

    cmake \
        -D MPI_DIR:FILEPATH=${MPI_PREFIX} \
        -D DISABLE_FORTRAN:BOOL=on \
        -D CMAKE_EXE_LINKER_FLAGS="-L${HDF5_PREFIX}/lib -lhdf5_hl -lhdf5 -L${MPI_PREFIX}/lib -lmpi -lmpi_cxx -L${ZLIB_PREFIX}/lib -lz " \
        -D NETCDF_DIR:PATH=${NETCDF_PREFIX} \
        -D CMAKE_INSTALL_PREFIX:PATH=${PREFIX} \
        ..

#       -D BUILD_TESTING:BOOL=off \

    if [ $? -ne 0 ]; then
        exit 
    fi
    make -j ${MAKE_NP}
    if [ $? -ne 0 ]; then
        exit 
    fi
    #make check
    if [ $? -ne 0 ]; then
        exit 
    fi
    make install
}

function build_exodus_make {
    export NETCDF_DIR=${NETCDF_PREFIX}

    EXODUS_DIR=${PREFIX}/exodusii/exodusii-${EXODUS_VERSION}
    rm -rf ${EXODUS_DIR}
    mkdir -p ${PREFIX}/exodusii
    tar ${TAR_FLAGS} ${DOWNLOAD_DIRECTORY}/exodusii-${EXODUS_VERSION}.tar.gz -C ${PREFIX}/exodusii

    cd ${EXODUS_DIR}
    cp -p Makefile.standalone Makefile.standalone.orig
    sed "s:NETCDF_INC = -I/usr/local/eng_sci/struct/x86_64/current64-gcc/inc:NETCDF_INC = -I$NETCDF_PREFIX/include:g" \
      < Makefile.standalone.orig > Makefile.standalone
    sed "s:NETCDF_LIB_DIR = /usr/local/eng_sci/struct/x86_64/current64-gcc/lib:NETCDF_LIB_DIR = $NETCDF_PREFIX/lib:g" \
      < Makefile.standalone > Makefile.standalone.1
    mv Makefile.standalone.1 Makefile.standalone

    make -j "$MAKE_NP" OS_TYPE=Linux BITS=64 CC=$MPICC CXX=$MPICXX -f Makefile.standalone
    cp -rp ./cbind/include ${EXODUS_PREFIX}/include
    cp -p *.a ${EXODUS_PREFIX}/lib
}

################################################################################
#
# metis
#
################################################################################
function build_metis {
    if [ ${METIS_PREFIX} == ${PREFIX} ]; then
        METIS_DIR=${PREFIX}/metis/metis-${METIS_VERSION}
        rm -rf ${METIS_DIR}
        mkdir -p ${PREFIX}/metis
        tar ${TAR_FLAGS} ${DOWNLOAD_DIRECTORY}/metis-${METIS_VERSION}.tar.gz -C ${PREFIX}/metis
        cd ${METIS_DIR}

        # Change CC to the mpi compiler
        mpicc_compiler=${MPICC}
        echo ">>>>${mpicc_compiler}"
        perl -w -i -p -e "s@^CC[\s]=.*@CC = ${mpicc_compiler}@" Makefile.in
    
        # no configuration...?
        make -j ${MAKE_NP}
        if [ $? -ne 0 ]; then
            exit 
        fi
        # need a manual install
        # copy the binary files:
        if [ ! -e ${METIS_PREFIX}/bin ]; then mkdir -p ${METIS_PREFIX}/bin ; fi
        cp pmetis kmetis oemetis onmetis partnmesh partdmesh mesh2nodal mesh2dual graphchk ${METIS_PREFIX}/bin
        # copy the libary file:
        if [ ! -e ${METIS_PREFIX}/lib ]; then mkdir -p ${METIS_PREFIX}/lib ; fi
        cp libmetis.a ${METIS_PREFIX}/lib
        # copy the include files
        if [ ! -e ${METIS_PREFIX}/include ]; then mkdir -p ${METIS_PREFIX}/include ; fi
        cp Lib/defs.h Lib/metis.h Lib/rename.h Lib/macros.h Lib/proto.h Lib/struct.h ${METIS_PREFIX}/include
    else
	echo Using METIS from ${METIS_PREFIX}
    fi

}

################################################################################
#
# Trilinos 
#
################################################################################
function build_trilinos {

    if [ ${TRILINOS_PREFIX} == $PREFIX ]; then

        TRILINOS_SRC=${PREFIX}/trilinos/trilinos-${TRILINOS_VERSION}-Source
        TRILINOS_BUILD=${PREFIX}/trilinos/trilinos-${TRILINOS_VERSION}-build
        TRILINOS_INSTALL=${PREFIX}/trilinos/trilinos-${TRILINOS_VERSION}-install 

        rm -rf ${TRILINOS_BUILD} ${TRILINOS_INSTALL} ${TRILINOS_SRC}
        mkdir -p ${PREFIX}/trilinos
        tar ${TAR_FLAGS} ${DOWNLOAD_DIRECTORY}/trilinos-${TRILINOS_VERSION}-Source.tar.gz -C ${PREFIX}/trilinos

        mkdir -p ${TRILINOS_BUILD}
        cd ${TRILINOS_BUILD}

	# Need this set so FindBoost.cmake works correctly
        export BOOST_ROOT=${BOOST_PREFIX}

        cmake \
                  -D CMAKE_BUILD_TYPE:STRING=DEBUG \
                  -D CMAKE_INSTALL_PREFIX:PATH=${TRILINOS_INSTALL} \
                  -D Trilinos_ENABLE_Fortran:BOOL=OFF \
                  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
                  -D Trilinos_ENABLE_STK:BOOL=ON \
                  -D Trilinos_ENABLE_TEUCHOS:BOOL=ON \
                  -D Trilinos_ENABLE_EPETRA:BOOL=ON \
                  -D Trilinos_ENABLE_BELOS:BOOL=ON \
                  -D Trilinos_ENABLE_NOX:BOOL=ON \
                  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON \
                  -D TPL_ENABLE_MPI:BOOL=ON \
		  -D TPL_ENABLE_ExodusII:BOOL=ON \
		  -D TPL_ExodusII_LIBRARIES:FILEPATH=${EXODUS_PREFIX}/lib \
		  -D TPL_ExodusII_INCLUDE_DIRS:FILEPATH=${EXODUS_PREFIX}/include \
                  -D Boost_INCLUDE_DIRS:FILEPATH=${BOOST_PREFIX}/include \
                  -D Boost_LIBRARY_DIRS:FILEPATH=${BOOST_PREFIX}/lib \
                  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF       \
                  -D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF   \
                  -D Trilinos_ENABLE_TESTS:BOOL=OFF \
                  -D CMAKE_C_COMPILER:FILEPATH=${MPICC} \
                  -D CMAKE_CXX_COMPILER:FILEPATH=${MPICXX} \
                  -D MPI_EXEC:FILEPATH=${MPIEXEC} \
                  ${TRILINOS_SRC}

        if [ $? -ne 0 ]; then
            exit 
        fi
        make -j${MAKE_NP}
        if [ $? -ne 0 ]; then
            exit 
        fi
        ctest -W 100
        if [ $? -ne 0 ]; then
            exit 
        fi
        make install

    else
	echo Using Trilinos found in ${TRILINOS_PREFIX}
    fi

}
################################################################################
#
# mstk
#
################################################################################
function build_mstk {
    if [ ${MSTK_PREFIX} == $PREFIX ]; then
        MSTK_DIR=${PREFIX}/mstk/mstk-${MSTK_VERSION}
        rm -rf ${MSTK_DIR}
        mkdir -p ${PREFIX}/mstk
        tar ${TAR_FLAGS} ${DOWNLOAD_DIRECTORY}/mstk-${MSTK_VERSION}.tgz -C ${PREFIX}/mstk
        #tar xvf ${DOWNLOAD_DIRECTORY}/mstk-${MSTK_VERSION}.tar -C ${PREFIX}/mstk
        cd ${MSTK_DIR}
   
	mkdir -p ${MSTK_DIR}/build
	cd ${MSTK_DIR}/build
        cmake \
    	    -D CMAKE_BUILD_TYPE:STRING=Release \
    	    -D ENABLE_PARALLEL=yes \
    	    -D ENABLE_ExodusII=yes \
	    -D HDF5_DIR:FILEPATH=${HDF5_PREFIX} \
    	    -D NetCDF_DIR:FILEPATH=${NETCDF_PREFIX} \
    	    -D ExodusII_DIR:FILEPATH=${EXODUS_PREFIX} \
    	    -D Metis_DIR:FILEPATH=${METIS_PREFIX} \
    	    -D METIS_LIB_DIR:FILEPATH=${METIS_PREFIX}/lib \
    	    -D Metis_INCLUDE_DIR:FILEPATH=${METIS_PREFIX}/include \
    	    -D ENABLE_Tests=no \
    	    -D INSTALL_DIR:FILEPATH=${MSTK_PREFIX} \
            -D INSTALL_ADD_VERSION=no \
	    ${MSTK_DIR}
    
            if [ $? -ne 0 ]; then
                exit 
            fi
            make -j ${MAKE_NP}
            if [ $? -ne 0 ]; then
                exit 
            fi
            make install

            # Copy the library 
	    if [ ! -e ${MSTK_PREFIX}/lib ]; then mkdir -p ${MSTK_PREFIX}/lib; fi
	    if [ -e libmstk.a ]; then cp libmstk.a ${MSTK_PREFIX}/lib; fi
	    if [ -e libmstk-d.a ]; then cp libmstk-d.a ${MSTK_PREFIX}/lib; fi
    else
	echo Using MSTK in ${MSTK_PREFIX}
    fi

}
################################################################################
#
# moab
#
################################################################################
function download_moab {
    # download moab from svn and tar it up
    mkdir -p ${PREFIX}/moab
    cd ${PREFIX}/moab
    NOW=`date "+%Y%m%d-%H%M%S"`
    svn co https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk moab-svn
    tar zcvf moab-svn-$NOW.tgz moab-svn
    mv moab-svn-$NOW.tgz $SOURCE
}

function build_moab {
    if [ $MOAB_PREFIX == $PREFIX ]; then
        MOAB_DIR=${MOAB_PREFIX}/moab/MOAB-${MOAB_VERSION}
	mkdir -p ${MOAB_PREFIX}/moab
	cd ${MOAB_PREFIX}/moab
	if [ $MOAB_VERSION == moab_svn ]; then # Go get the latest copy
            NOW=`date "+%Y%m%d-%H%M%S"`
            svn co https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk MOAB-${MOAB_VERSION}
	    cd MOAB-${MOAB_VERSION}
	    autoreconf -fi
	else 
	    tarfile=${DOWNLOAD_DIRECTORY}/MOAB-${MOAB_VERSION}.tar.gz
	    tar ${TAR_FLAGS} $tarfile 
	    mv MOAB MOAB-${MOAB_VERSION}
	fi
        cd ${MOAB_DIR}
    
        # need to uncomment this if compiling tests
        #perl -w -i -p -e "s@malloc\.h@stdlib.h@" ${MOAB_DIR}/itaps/imesh/testc_cbind.c

        export LDFLAGS="-L${MPI_PREFIX}/lib -lmpi -L${HDF5_PREFIX}/lib -lhdf5_hl -lhdf5 -L${ZLIB_PREFIX}/lib -lz"

        ./configure --prefix=${PREFIX} \
                    --disable-fortran \
                    --with-mpi=${MPI_PREFIX} \
                    --with-hdf5=${HDF5_PREFIX} \
                    --with-netcdf=${NETCDF_PREFIX} 

         if [ $? -ne 0 ]; then
            exit 
         fi
         make -j ${MAKE_NP}
         if [ $? -ne 0 ]; then
            exit 
         fi
         # moab tests don't compile...
         #make check
         #if [ $? -ne 0 ]; then
         #   exit 
         #fi
         make install
     else
	 echo Using MOAB found in ${MOAB_PREFIX}
     fi

}

################################################################################
#
# write a script for configuring and building amanzi
#
################################################################################
function generate_amanzi_build_script {
    cd ${SCRIPT_DIR}
    if [ -e ${AMANZI_BUILD_SCRIPT} ]; then mv ${AMANZI_BUILD_SCRIPT} ${AMANZI_BUILD_SCRIPT}.save; fi
    cat > ${AMANZI_BUILD_SCRIPT} <<EOF
#!/bin/bash

AMANZI_WORK_DIR=${WORK_DIR}
AMANZI_CONFIG=0
AMANZI_DIR=
AMANZI_MAKE=0
AMANZI_CLOBBER=0
AMANZI_TEST=0
PLATFORM=${PLATFORM}

# Add module command if available
if [[ \$MODULESHOME ]]; then
    . \$MODULESHOME/init/bash
fi

# If runnnig on NERSC machines (Franklin, Hopper, etc.), use Amanzi module.
if [[ \$NERSC_HOST ]]; then
    module use /project/projectdirs/m1012/modulefiles/\$NERSC_HOST
    
    if [[ \$PE_ENV ]]; then
	compiler_loaded=`echo \$PE_ENV | tr "[A-Z]" "[a-z]"`
	module unload PrgEnv-\${compiler_loaded}
    fi

    module load AmanziEnv-gnu

fi

function determine_amanzi_dir() {
    AMANZI_PATH=
    case \$1 in
	/*) AMANZI_DIR=\$1;;
	*) AMANZI_DIR=\${AMANZI_WORK_DIR}/\$1;;
    esac
    #echo \$AMANZI_DIR
}

# Process the command line arguments
while getopts "abcd:mt" flag
do
  case \$flag in
    a) AMANZI_CONFIG=1; AMANZI_MAKE=1; AMANZI_TEST=1;;
    b) AMANZI_CLOBBER=1;;
    c) AMANZI_CONFIG=1;;
    d) determine_amanzi_dir \${OPTARG};;
    m) AMANZI_MAKE=1;;
    t) AMANZI_TEST=1;;
  esac
done

echo AMANZI_CONFIG=\$AMANZI_CONFIG
echo AMANZI_DIR=\$AMANZI_DIR
echo AMANZI_MAKE=\$AMANZI_MAKE
echo AMANZI_TEST=\$AMANZI_TEST

if [ \${PLATFORM} == Fedora ]; then
    LD_RUN_PATH=\$LD_RUN_PATH:${PREFIX}/lib
    LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:${PREFIX}/lib
    export LD_RUN_PATH LD_LIBRARY_PATH
fi

if [ \$AMANZI_CLOBBER -eq 1 ]; then
    rm -rf \${AMANZI_DIR}/build
    cd \${AMANZI_DIR}/src
    if [ -f ./Makefile ]; then
        make clean
    fi
    rm -rf \\
        CMakeCache.txt \\
        *.cmake */*.cmake */*/*.cmake \\
        CMakeFiles */CMakeFiles */*/CMakeFiles \\
        Makefile */Makefile */*/Makefile \\
        Testing \\
        *~ */*~ */*/*~ \\
        *.rej */*.rej */*/*.rej
fi


if [ \$AMANZI_CONFIG -eq 1 ]; then
    # for out of source builds
    #rm -rf \${AMANZI_DIR}/build
    mkdir -p \${AMANZI_DIR}/build
    cd \${AMANZI_DIR}/build
    # for in source builds
#    cd \${AMANZI_DIR}/src

    export CXX=${MPICXX}
    export CC=${MPICC}
    export BOOST_ROOT=${BOOST_PREFIX}

    cmake \\
        -D ENABLE_Config_Report:BOOL=ON \\
        -D MPI_DIR:FILEPATH=${MPI_PREFIX} \\
        -D MPI_EXEC:FILEPATH=${MPIEXEC} \\
        -D MPI_EXEC_NUMPROCS_FLAG:STRING=${MPIEXEC_NUMPROCS_FLAG} \\
        -D MPI_EXEC_MAX_NUMPROCS:STRING=${MPIEXEC_MAX_NUMPROCS} \\
        -D ENABLE_TESTS:BOOL=ON \\
        -D UnitTest_DIR:FILEPATH=${UNITTEST_PREFIX} \\
        -D HDF5_DIR:FILEPATH=${HDF5_PREFIX} \\
		-D ENABLE_ASCEMIO:BOOL=ON \\
		-D ASCEMIO_DIR=${ASCEMIO_PREFIX} \\
        -D NetCDF_DIR:FILEPATH=${NETCDF_PREFIX} \\
        -D ExodusII_DIR:FILEPATH=${EXODUS_PREFIX} \\
        -D ENABLE_MESH_MOAB:BOOL=ON \\
        -D MOAB_DIR:FILEPATH=${MOAB_PREFIX} \\
        -D ENABLE_MESH_MSTK:BOOL=ON \\
        -D MSTK_DIR:FILEPATH=${MSTK_PREFIX} \\
        -D METIS_DIR:FILEPATH=${METIS_PREFIX} \\
        -D ENABLE_STK_Mesh:BOOL=ON \\
        -D Trilinos_DIR:FILEPATH=${TRILINOS_PREFIX}/trilinos-${TRILINOS_VERSION}-install \\
	-D ENABLE_Unstructured:Bool=ON \\
	-D ENABLE_Structured:Bool=ON \\
        ..
#        ../src
# for out of source: ../source    for in source: .
    if [ \$? -ne 0 ]; then
        exit 1
    fi
fi

if [ \$AMANZI_MAKE -eq 1 ]; then
    cd \${AMANZI_DIR}/build
#    cd \${AMANZI_DIR}/src
    make -j ${MAKE_NP}
    if [ \$? -ne 0 ]; then
        exit 1
    fi
fi

if [ \$AMANZI_TEST -eq 1 ]; then
    cd \${AMANZI_DIR}/build
#    cd \${AMANZI_DIR}/src
    ctest --timeout 60 --output-on-failure
    if [ \$? -ne 0 ]; then
        exit 1
    fi
fi

EOF
    chmod 700 amanzi-build.sh
}



################################################################################
#
# MAIN
#
#################################################################################

# Flag to print the help message
PRINT_HELP=0

# Found errors on the command line
OPTS_OK=1

# Process command line
while getopts "hd:aAbBegHkmnstuwz" flag
do
  case $flag in
      h) PRINT_HELP=1;;
      d) DOWNLOAD_DIRECTORY=${OPTARG};;
      a) BUILD_BOOST=1; BUILD_EXODUS=1; BUILD_HDF5=1; BUILD_MSTK=1; BUILD_MOAB=1; BUILD_NETCDF=1; BUILD_METIS=1; BUILD_TRILINOS=1; BUILD_UNITTEST=1; BUILD_CURL=1; BUILD_ZLIB=1;;
	  A) BUILD_ASCEMIO=1;;
      b) BUILD_BOOST=1;;
      B) BUILD_AMANZI_SCRIPT=1;;
      e) BUILD_EXODUS=1;;
      H) BUILD_HDF5=1;;
      k) BUILD_MSTK=1;;
      m) BUILD_MOAB=1;;
      n) BUILD_NETCDF=1;;
      s) BUILD_METIS=1;;
      t) BUILD_TRILINOS=1;;
      u) BUILD_UNITTEST=1;;
      w) BUILD_CURL=1;;
      z) BUILD_ZLIB=1;;
      \?)
          echo " Invalid option $OPTARG"
	  OPTS_OK=0
	  ;;
      :) 
          echo "$OPTARG requires an argument"
	  OPTS_OK=0
	  ;;
  esac
done

if [ $PRINT_HELP == 1 ]; then
    help_message
fi

if [ ! $OPTS_OK ]; then
    help_message
fi

print_build_status
check_platform
check_mpi
check_download
check_prefix

#
# build non-mpi dependent libraries and mpi using system compilers
#
echo "Building non-MPI dependent libraries"
if [ $BUILD_CURL -eq 1 ]; then
    build_curl
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_ZLIB -eq 1 ]; then
    build_zlib
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_UNITTEST -eq 1 ]; then
    build_unittest
    cd ${SCRIPT_DIR}
fi

#
#  build mpi dependent libraries using mpi compilers
#
echo "Building MPI dependent libraries"

export CXX=${MPICXX}
export CC=${MPICC}

if [ $BUILD_BOOST -eq 1 ]; then
    build_boost
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_HDF5 -eq 1 ]; then
    build_hdf5
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_ASCEMIO -eq 1 ]; then
    build_ascemio
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_NETCDF -eq 1 ]; then
    build_netcdf
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_EXODUS -eq 1 ]; then
    build_exodus_cmake
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_METIS -eq 1 ]; then
    build_metis
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_TRILINOS -eq 1 ]; then
    build_trilinos
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_MSTK -eq 1 ]; then
    build_mstk
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_MOAB -eq 1 ]; then
    build_moab
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_AMANZI_SCRIPT -eq 1 ]; then
    generate_amanzi_build_script
    echo "Writing ${AMANZI_BUILD_SCRIPT}"
fi
