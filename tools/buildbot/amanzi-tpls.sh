#!/bin/bash

################################################################################
#
# script to help automate building all the dependancies for amanzi
#
################################################################################
#
# Author: Ben Andre <bandre@lbl.gov>
#
################################################################################
#
# Instructions:
#
#   - Feel free to copy this file and modify it as needed to get all
#     the TPLs built on your system. Please do not commit changes to
#     this file unless you have tested it on all supported
#     systems. Small changes for one system can break things on
#     another system in a way that is difficult to debug.
#
#   - select one of the example configuration files,
#     amanzi-tpls.XXXXX, that is most appropriate for your system.
#
#   - Use your package manager to install: wget, boost, curl, unittest++, etc
#
#   - Optional: Use the package manager to install standard versions of: netcdf, hdf5-18
#
#   - Download the source for the packages you want installed by this script.
#     Trilinos and moab must be downloaded manually (download sites are listed
#     below). The remaining TPLs can be automatically downloaded with the "-d"
#     option to amanzi-tpls.sh.
#
#   - Place the package archives into a single directory
#
#   - copy your configuration file to ${HOME}/.amanzi-tpls and open it in an editor.
#
#   - Edit the PREFIX variable to be the directory you want to install into.
#
#   - Edit the WORK_DIR variable to be the directory containing your amanzi repositories
#
#   - Edit the SOURCE variable to be the directory where you saved the source.
#
#   - Edit the PACKAGE_NAME_PREFIX variables to point to the location for each
#     package. If you are using a system package, it should point to the system
#     location. If you are installing the package from this script, it should
#     point to the PREFIX directory. For example:
#
#     MPI_PREFIX=/usr
#     NetCDF_PREFIX=${PREFIX}
#
#   - Edit the PARALLEL_NP variable to be the number of parallel processes you want
#     to use to for building (make -j 4) and running MPI (mpiexec -np 4)
#
#   - To download, configure and build all packages, run this script as:
#
#     $ ./amanzi-libs.sh -d -a
#
#   - Individual packages can be built by using the appropriate
#     commandline options. See below.
#
#   - Wait.... The script should stop if any problems arise while
#     configuring, building or installing a package.
#
#   - this script creates another script called "amanzi-build.sh". This script
#     is used to configure, build and test amanzi. The -d flag is used to specify
#     which directory to work in. For example, "-d amanzi-debug" will look for an
#     amanzi repository within ${WORK_DIR}/amanzi-debug. Additional flags
#     control the build/configure/test/cleanup process.
#     Key flags are:
#          -a   : performs configure/build/test
#          -b   : "clobbers" the repository, tries to remove all generated and
#                  CMake files, forces rebuild of everything
#          -c   : run cmake in the repository
#          -d   : repository name
#          -m   : runs make in the repository
#          -t   : runs ctest in the repository
#
#     Examples:
#
#       - Obtain a fresh repository and build everything:
#         $ ./amanzi-build -d amanzi-debug -a -b
#
#       - Configure/make/test:
#         $ ./amanzi-build -d amanzi-debug -a
#
#       - recompile modified files and run the tests:
#         $ ./amanzi-build -d amanzi-debug -m -t
#
#   - Usage Notes:
#
#     - This does out of source builds by default. If you modify an input
#       file, you need to rerun the configure process to copy the input file
#       into the build directory.
#
#     - You can change paths in the configuration file and regenerate
#       just the amanzi-build.sh script by running amanzi-tpls.sh
#       with no options:
#       $ ./amanzi-tpls.sh
#
################################################################################
#
# PREFIX : where the packages are installed, e.g. ${HOME}/local
# SOURCE : where the packages archives are saved, e.g. ${HOME}/source
# MPI_PREFIX : where mpi lives, e.g. /usr, ${HOME}/local
#
################################################################################
#
# Notes:
#
#  - use the package manager to install everything it can or install
#    everything manually
#
#  - default mpi from OSX does not have fortran versions of the
#    compilers. if you want/need fortran support, mpi must be
#    installed through the package manager or build manually.
#
#  - the built in compilers on OSX 10.6 are gcc 4.2.1. mac ports
#    version of gcc 4.2 won't build on 10.6, so you must use a more
#    recent version of gcc. More recent version of gcc are not binary
#    compatible with gcc 4.2, so if you try linking to gcc 4.2
#    compiled libraries with gcc 4.4 compiled libraries, link time
#    errors will occur. The only way I've found to get around this is
#    install gcc 4.4-4.5 with mac ports, then compile all additional
#    libraries from scratch with this set of compilers. Installing mpi
#    with macports puts the executables at /opt/local/bin/openmpicxx
#    instead of /opt/local/bin/mpicxx. It is easier to maintain
#    compatability in the script if we just compile mpi here as well.
#
#  - build files from an existing build by this script are removed.
#    if a package is already installed in $PREFIX, the installed files
#    are not uninstalled!
#
#  - The rm -rf on source/build directories is intentional to start
#    with a fresh copy of the source every time. It helps with
#    automation and debugging, but increases build time.
#
#  - exodusii requires a hand modified version of netcdf to use
#    large meshes. using the package manager to install netcdf will
#    most likeley not have this modification.
#
################################################################################
#
# TODO:
#
#  - add more error checking to verify that configure/cmake/make did not fail....
#
#  - integrate the do_configure stripts from the amanzi repo...?
#
#  - flag to control in source/out of source builds
#
#  - HDF5 test is very long and locks up the machine....
#
#  - check the platform, look for a package manager and run the package
#    manager automatically?
#
#################################################################################
#
# Automatic downloads:
#
# The supported version of the TPLs can be automatically
# downloaded if you have wget on your system by calling the
# download_archives function (-d). They will be saved into ${SOURCE}.
#
#################################################################################
#
# TPL project home pages:
#
# openmpi:
# http://www.open-mpi.org
#
# boost:
# http://boost.org
#
# unittest++:
# http://sourceforge.net/projects/unittest-cpp
#
# zlib:
# http://zlib.net
#
# curl:
# http://curl.haxx.se
#
# hdf5:
# http://www.hdfgroup.org
#
# netcdf:
# http://www.unidata.ucar.edu
#
# exodusii:
# http://sourceforge.net/projects/exodusii
#
# metis:
# http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-4.0.3.tar.gz
#
# mstk:
# https://software.lanl.gov/MeshTools/trac
#
# ccse:
# https://ccse.lbl.gov/Software/tarfiles/ccse-0.1.7.tar.gz
#
# trilinos:
# http://trilinos.sandia.gov
#
# ascemio:
#
#
#
# moab:
#
#
#
#
################################################################################


################################################################################
#
# Very little below this point should need to be changed....
#
#################################################################################

SCRIPT_DIR=$PWD

TAR_FLAGS=xf

OPENMPI_VERSION=1.5
BOOST_VERSION=1_46_1
UNITTEST_VERSION=1.4
ZLIB_VERSION=1.2.5
CURL_VERSION=7.21.2
HDF5_VERSION=1.8.8
NETCDF_VERSION=4.1.1
EXODUS_VERSION=4.98
MOAB_VERSION=r4276
METIS_VERSION=4.0.3
MSTK_VERSION=1.86rc2
TRILINOS_VERSION=10.6.4
CCSE_VERSION=1.0.9
ASCEMIO_VERSION=2.1

AMANZI_TPL_ARCHIVES=https://software.lanl.gov/ascem/tpls
WGET_FLAGS=--no-check-certificate

################################################################################
#
# check that we are on a tested platform
#
################################################################################
function check_platform {
    PLATFORM=
    OS=`uname -s`
    if [ $OS == Darwin ]; then
	PLATFORM=Darwin
    elif [ $OS == Linux ]; then
	if [ -e /etc/lsb-release ]; then
	    # probably ubuntu.... may need to grep the file?
	    PLATFORM=Ubuntu
	elif [ -e /etc/fedora-release ]; then
	    PLATFORM=Fedora
	elif [ -e /etc/redhat-release ]; then
	    # centos and redhat are suppose to be the same....
	    PLATFORM=RHEL
	fi
    fi

    if [ $PLATFORM ]; then
	echo "Running on a $PLATFORM system..."
    else
	echo "You are running on an untested operating system ($OS, $PLATFORM). Good luck."
    fi
}

################################################################################
#
# download publically available archives with wget
#
################################################################################
function download_archives {
    mkdir -p ${SOURCE}
    cd ${SOURCE}
    if [ ${MPI_PREFIX} == ${PREFIX} -a ! -f ${SOURCE}/openmpi-${OPENMPI_VERSION}.tar.gz ]; then
	URL=${AMANZI_TPL_ARCHIVES}/openmpi-${OPENMPI_VERSION}.tar.gz
	wget ${WGET_FLAGS} $URL
    fi

    if [ ${BOOST_PREFIX} == ${PREFIX} -a ! -f ${SOURCE}/boost_${BOOST_VERSION}.tar.bz2 ]; then
	URL=${AMANZI_TPL_ARCHIVES}/boost_${BOOST_VERSION}.tar.bz2
	wget ${WGET_FLAGS} $URL
    fi

    if [ ${UNITTEST_PREFIX} == ${PREFIX} -a ! -f ${SOURCE}/unittest-cpp-${UNITTEST_VERSION}.zip ]; then
	URL=${AMANZI_TPL_ARCHIVES}/unittest-cpp-${UNITTEST_VERSION}.zip
	wget ${WGET_FLAGS} $URL
    fi

    if [ ${ZLIB_PREFIX} == ${PREFIX} -a ! -f ${SOURCE}/zlib-${ZLIB_VERSION}.tar.gz ]; then
	URL=${AMANZI_TPL_ARCHIVES}/zlib-${ZLIB_VERSION}.tar.gz
	wget ${WGET_FLAGS} $URL
    fi

    if [  ${CURL_PREFIX} == ${PREFIX} -a ! -f ${SOURCE}/curl-${CURL_VERSION}.tar.gz ]; then
	URL=${AMANZI_TPL_ARCHIVES}/curl-${CURL_VERSION}.tar.gz
	wget ${WGET_FLAGS} $URL
    fi

    if [  ${HDF5_PREFIX} == ${PREFIX} -a ! -f ${SOURCE}/hdf5-${HDF5_VERSION}.tar.gz ]; then
	URL=${AMANZI_TPL_ARCHIVES}/hdf5-${HDF5_VERSION}.tar.gz
	wget ${WGET_FLAGS} $URL
    fi

    if [  ${NETCDF_PREFIX} == ${PREFIX} -a ! -f ${SOURCE}/netcdf-${NETCDF_VERSION}.tar.gz ]; then
	URL=${AMANZI_TPL_ARCHIVES}/netcdf-${NETCDF_VERSION}.tar.gz
	wget ${WGET_FLAGS} $URL
    fi

    if [  ! -f ${SOURCE}/exodusii-${EXODUS_VERSION}.tar.gz ]; then
	URL=${AMANZI_TPL_ARCHIVES}/exodusii-${EXODUS_VERSION}.tar.gz
	wget ${WGET_FLAGS} $URL
    fi

    if [  ! -f ${SOURCE}/metis-${METIS_VERSION}.tar.gz ]; then
	URL=${AMANZI_TPL_ARCHIVES}/metis-${METIS_VERSION}.tar.gz
	wget ${WGET_FLAGS} $URL
    fi

    if [  ! -f ${SOURCE}/mstk-${MSTK_VERSION}.tgz ]; then
	URL=${AMANZI_TPL_ARCHIVES}/mstk-${MSTK_VERSION}.tgz
	wget ${WGET_FLAGS} $URL
    fi

    if [  ! -f ${SOURCE}/ccse-${CCSE_VERSION}.tar.gz ]; then
	URL=${AMANZI_TPL_ARCHIVES}/ccse-${CCSE_VERSION}.tar.gz
	wget ${WGET_FLAGS} $URL
    fi

    if [  ! -f ${SOURCE}/ascem-io-${ASCEMIO_VERSION}.tar.gz ]; then
	URL=${AMANZI_TPL_ARCHIVES}/ascem-io-${ASCEMIO_VERSION}.tar.gz
	wget ${WGET_FLAGS} $URL
    fi

    if [  ! -f ${SOURCE}/MOAB-${MOAB_VERSION}.tar.gz ]; then
	URL=${AMANZI_TPL_ARCHIVES}/MOAB-${MOAB_VERSION}.tar.gz
	wget ${WGET_FLAGS} $URL
    fi

    if [  ! -f ${SOURCE}/trilinos-${TRILINOS_VERSION}-Source.tar.gz ]; then
	URL=${AMANZI_TPL_ARCHIVES}/trilinos-${TRILINOS_VERSION}-Source.tar.gz
	wget ${WGET_FLAGS} $URL
    fi

    cd ${SCRIPT_DIR}
}

################################################################################
#
# package managers....
#
################################################################################

# NOTE: in macports, netcdf depends on hdf5-18, not hdf5.
function install_macports {
    sudo port install boost
    sudo port install unittest-cpp
#    sudo port install hdf5-18 +openmpi -universal
#    sudo port install netcdf +openmpi -universal
}

function install_apt_get {
    sudo apt-get install XXXXXXX
}

################################################################################
#
# mpi
#
################################################################################
function install_mpi {
    if [ ${MPI_PREFIX} == ${PREFIX} ]; then
	MPI_DIR=${PREFIX}/mpi/openmpi-${OPENMPI_VERSION}
	rm -rf ${MPI_DIR}
	mkdir -p ${PREFIX}/mpi
	tar ${TAR_FLAGS} ${SOURCE}/openmpi-${OPENMPI_VERSION}.tar.gz -C ${PREFIX}/mpi
	cd ${MPI_DIR}
	./configure --prefix=${PREFIX} \
	    --enable-static

	if [ $? -ne 0 ]; then
	    exit
	fi
	make -j${PARALLEL_NP} all
	if [ $? -ne 0 ]; then
	    exit
	fi
	make install
	if [ $? -ne 0 ]; then
	    exit
	fi
    else
	echo Using mpi from ${MPI_PREFIX}
    fi
}

################################################################################
#
# boost
#
# only build the libraries that amanzi explicitly depends on
#
################################################################################
function install_boost {
    if [ ${BOOST_PREFIX} == ${PREFIX} ]; then
	BOOST_DIR=${PREFIX}/boost/boost_${BOOST_VERSION}
	rm -rf ${BOOST_DIR}
	mkdir -p ${PREFIX}/boost
	tar ${TAR_FLAGS} ${SOURCE}/boost_${BOOST_VERSION}.tar.bz2 -C ${PREFIX}/boost
	cd ${BOOST_DIR}

	./bootstrap.sh --prefix=${PREFIX} \
	    --with-libraries=system,filesystem,program_options,regex
	if [ $? -ne 0 ]; then
	    exit
	fi
	./bjam
	if [ $? -ne 0 ]; then
	    exit
	fi
	./bjam install
	if [ $? -ne 0 ]; then
	    exit
	fi
    else
	echo Using Boost from ${BOOST_PREFIX}
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
function install_zlib {
    if [ ${ZLIB_PREFIX} == ${PREFIX} ]; then
	ZLIB_DIR=${PREFIX}/zlib/zlib-${ZLIB_VERSION}
	rm -rf ${ZLIB_DIR}
	mkdir -p ${PREFIX}/zlib
	tar ${TAR_FLAGS} ${SOURCE}/zlib-${ZLIB_VERSION}.tar.gz -C ${PREFIX}/zlib
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
# curl
#
################################################################################
function install_curl {
    if [ ${CURL_PREFIX} == ${PREFIX} ]; then
	CURL_DIR=${PREFIX}/curl/curl-${CURL_VERSION}
	rm -rf ${CURL_DIR}
	mkdir -p ${PREFIX}/curl
	tar ${TAR_FLAGS} ${SOURCE}/curl-${CURL_VERSION}.tar.gz -C ${PREFIX}/curl
	cd ${CURL_DIR}
	./configure --prefix=${PREFIX} --enable-static
	if [ $? -ne 0 ]; then
	    exit
	fi
	make -j${PARALLEL_NP} all
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
# unittest++
#
################################################################################
function install_unittest {
    if [ ${UNITTEST_PREFIX} == ${PREFIX} ]; then
	UNITTEST_DIR=${PREFIX}/unittest/UnitTest++
	rm -rf ${UNITTEST_DIR}
	mkdir -p ${PREFIX}/unittest
	unzip ${SOURCE}/unittest-cpp-${UNITTEST_VERSION}.zip -d ${PREFIX}/unittest

	cd ${UNITTEST_DIR}
	make -j${PARALLEL_NP} all
	if [ $? -ne 0 ]; then
	    exit
	fi
	# ugh... manual install
	echo Copying UnitTest++ files to ${UNITTEST_PREFIX}
	cp libUnitTest++.a ${UNITTEST_PREFIX}/lib
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
# hdf5
#
# parallel is not compatible with c++
# (http://www.hdfgroup.org/hdf5-quest.html#p5thread) so which is necessary
# for amanzi?
#
################################################################################
function install_hdf5 {
    if [ ${HDF5_PREFIX} == ${PREFIX} ]; then
	HDF5_DIR=${PREFIX}/hdf5/hdf5-${HDF5_VERSION}
	rm -rf ${HDF5_DIR}
	mkdir -p ${PREFIX}/hdf5
	tar ${TAR_FLAGS} ${SOURCE}/hdf5-${HDF5_VERSION}.tar.gz -C ${PREFIX}/hdf5

	cd ${HDF5_DIR}

	./configure --prefix=${PREFIX} \
	    --disable-fortran \
	    --enable-production \
	    --enable-largefile \
	    --enable-parallel

	if [ $? -ne 0 ]; then
	    exit
	fi
	make -j ${PARALLEL_NP}
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
# netcdf
#
################################################################################
function install_netcdf {
    if [ ${NETCDF_PREFIX} == ${PREFIX} ]; then
	NETCDF_DIR=${PREFIX}/netcdf/netcdf-${NETCDF_VERSION}
	rm -rf ${NETCDF_DIR}
	mkdir -p ${PREFIX}/netcdf
	tar ${TAR_FLAGS} ${SOURCE}/netcdf-${NETCDF_VERSION}.tar.gz -C ${PREFIX}/netcdf

	cd ${NETCDF_DIR}

	perl -w -i -p -e 's@#define NC_MAX_DIMS[\s]+[\d]+@#define NC_MAX_DIMS 65536@' libsrc/netcdf.h libsrc4/netcdf.h libsrc4/netcdf_base.h
	perl -w -i -p -e 's@#define NC_MAX_VARS[\s]+[\d]+@#define NC_MAX_VARS 524288@' libsrc/netcdf.h libsrc4/netcdf.h libsrc4/netcdf_base.h
	perl -w -i -p -e 's@#define NC_MAX_VAR_DIMS[\s]+NC_MAX_DIMS@#define NC_MAX_VAR_DIMS 8@' libsrc/netcdf.h libsrc4/netcdf.h libsrc4/netcdf_base.h

#        grep -e NC_MAX_DIMS libsrc4/netcdf.h
#        grep -e NC_MAX_VARS libsrc4/netcdf.h
#        grep -e NC_MAX_VAR_DIMS libsrc4/netcdf.h

	USE_NETCDF4='--disable-netcdf-4 --disable-cxx-4'
	if [ $PLATFORM == Ubuntu ]; then
	    USE_NETCDF4='--enable-netcdf-4 --enable-cxx-4'
	fi

	./configure --prefix=${PREFIX} \
	    --disable-fortran \
	    --disable-f90 \
	    --disable-f77 \
	    --disable-fortran-compiler-check \
	    ${USE_NETCDF4} \
	    --disable-dap \
	    --with-mpi=${MPI_PREFIX} \
	    --with-hdf5=${HDF5_PREFIX} \
	    --with-zlib=${ZLIB_PREFIX}
	if [ $? -ne 0 ]; then
	    exit
	fi
	make -j ${PARALLEL_NP}
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
function install_exodus {
    # download source from: http://sourceforge.net/projects/exodusii/
    export NETCDF_DIR=${NETCDF_PREFIX}

    EXODUS_DIR=${PREFIX}/exodusii/exodusii-${EXODUS_VERSION}
    rm -rf ${EXODUS_DIR}
    mkdir -p ${PREFIX}/exodusii
    tar ${TAR_FLAGS} ${SOURCE}/exodusii-${EXODUS_VERSION}.tar.gz -C ${PREFIX}/exodusii
    mkdir -p ${EXODUS_DIR}/build


    cd ${EXODUS_DIR}

    perl -w -i -p -e 's@TARGET_LINK_LIBRARIES\(exoIIv2c \${NETCDF_LIBRARY}/libnetcdf\.a\)@TARGET_LINK_LIBRARIES\(exoIIv2c \${NETCDF_DIR}/lib/libnetcdf\.a \${HDF5_PREFIX}/lib/libhdf5_hl\.a \${HDF5_PREFIX}/lib/libhdf5\.a\)@' cbind/CMakeLists.txt

    cd ${EXODUS_DIR}/build

    cmake \
	-D NETCDF_DIR:FILEPATH=${NETCDF_DIR} \
	-D HDF5_PREFIX:PATH=${HDF5_PREFIX} \
	-D CMAKE_EXE_LINKER_FLAGS="-L${HDF5_PREFIX}/lib -lhdf5_hl -lhdf5 -L${MPI_PREFIX}/lib -lmpi -lmpi_cxx -L${ZLIB_PREFIX}/lib -lz " \
	-D CMAKE_INSTALL_PREFIX:PATH=${PREFIX} \
	..

#       -D BUILD_TESTING:BOOL=off \

    if [ $? -ne 0 ]; then
	exit
    fi
    make -j ${PARALLEL_NP}
    if [ $? -ne 0 ]; then
	exit
    fi
    #make check
    if [ $? -ne 0 ]; then
	exit
    fi
    make install
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

function install_moab {
    MOAB_DIR=${PREFIX}/moab/MOAB
    rm -rf ${MOAB_DIR}
    mkdir -p ${PREFIX}/moab
    tar ${TAR_FLAGS} ${SOURCE}/MOAB-${MOAB_VERSION}.tar.gz -C ${PREFIX}/moab
    cd ${MOAB_DIR}

    # need to uncomment this if compiling tests
    #perl -w -i -p -e "s@malloc\.h@stdlib.h@" ${MOAB_DIR}/itaps/imesh/testc_cbind.c

    #autoreconf -fi
    export LDFLAGS="-L${MPI_PREFIX}/lib -lmpi -L${HDF5_PREFIX}/lib -lhdf5_hl -lhdf5 -L${ZLIB_PREFIX}/lib -lz"

    ./configure --prefix=${PREFIX} \
	--disable-fortran \
	--with-mpi=${MPI_PREFIX} \
	--with-hdf5=${HDF5_PREFIX} \
	--with-netcdf=${NETCDF_PREFIX}

    if [ $? -ne 0 ]; then
	exit
    fi
    make -j ${PARALLEL_NP}
    if [ $? -ne 0 ]; then
	exit
    fi
    # moab tests don't compile...
    #make check
    #if [ $? -ne 0 ]; then
    #   exit
    #fi
    make install

}

################################################################################
#
# metis
#
################################################################################
function install_metis {
    METIS_DIR=${PREFIX}/metis/metis-${METIS_VERSION}
    rm -rf ${METIS_DIR}
    mkdir -p ${PREFIX}/metis
    tar ${TAR_FLAGS} ${SOURCE}/metis-${METIS_VERSION}.tar.gz -C ${PREFIX}/metis
    cd ${METIS_DIR}

    # Change CC to the mpi compiler
    perl -w -i -p -e "s@^CC[\s]=.*@CC = ${CC}@" Makefile.in

    # no configuration...?
    make -j ${PARALLEL_NP}
    if [ $? -ne 0 ]; then
	exit
    fi
    # need a manual install
    # copy the binary files:
    cp pmetis kmetis oemetis onmetis partnmesh partdmesh mesh2nodal mesh2dual graphchk ${METIS_PREFIX}/bin
    # copy the libary file:
    cp libmetis.a ${METIS_PREFIX}/lib
    # copy the include files
    cp Lib/defs.h Lib/metis.h Lib/rename.h Lib/macros.h Lib/proto.h Lib/struct.h ${METIS_PREFIX}/include
}

################################################################################
#
# ascemio
#
################################################################################
function install_ascemio {
    ASCEMIO_DIR=${PREFIX}/ascem-io/ascem-io-${ASCEMIO_VERSION}
    rm -rf ${ASCEMIO_DIR}
    mkdir -p ${PREFIX}/ascem-io
    tar ${TAR_FLAGS} ${SOURCE}/ascem-io-${ASCEMIO_VERSION}.tar.gz -C ${PREFIX}/ascem-io
    cd ${ASCEMIO_DIR}

    # no configuration, just make
    cd ${ASCEMIO_DIR}/src
    make CC=${CC} \
	HDF5_INCLUDE_DIR=${HDF5_PREFIX}/include

    if [ $? -ne 0 ]; then
	exit
    fi
    # testing...?

    # install
    make ASCEMIO_INSTALL_DIR=${PREFIX} install
    if [ $? -ne 0 ]; then
	exit
    fi
}

################################################################################
#
# mstk
#
################################################################################
function install_mstk {
    MSTK_DIR=${PREFIX}/mstk/mstk-${MSTK_VERSION}
    rm -rf ${PREFIX}/mstk
    mkdir -p ${PREFIX}/mstk
    tar ${TAR_FLAGS} ${SOURCE}/mstk-${MSTK_VERSION}.tgz -C ${PREFIX}/mstk
    cd ${MSTK_DIR}

    # for some reason installing examples dies on centos with cmake 2.8.4
    perl -w -i -p -e 's/(install\(TARGETS serial_example)/#$1/' ${MSTK_DIR}/example/CMakeLists.txt
    perl -w -i -p -e 's/(install\(TARGETS parallel_example)/#$1/' ${MSTK_DIR}/example/CMakeLists.txt

    cmake \
	-Wno-dev \
	-D CMAKE_BUILD_TYPE:STRING=Release \
	-D ENABLE_PARALLEL=yes \
	-D ENABLE_ExodusII=yes \
        -D ENABLE_ZOLTAN:BOOL=TRUE \
        -D ZOLTAN_DIR:PATH=${TRILINOS_INSTALL} \
	-D HDF5_DIR:FILEPATH=${HDF5_PREFIX} \
	-D NetCDF_DIR:FILEPATH=${NETCDF_PREFIX} \
	-D ExodusII_DIR:FILEPATH=${PREFIX} \
	-D Metis_DIR:FILEPATH=${METIS_PREFIX} \
	-D METIS_LIB_DIR:FILEPATH=${METIS_PREFIX}/lib \
	-D Metis_INCLUDE_DIR:FILEPATH=${METIS_PREFIX}/include \
	-D ENABLE_Tests=no \
	-D INSTALL_DIR:FILEPATH=${MSTK_PREFIX} \
	-D INSTALL_ADD_VERSION=no \
	./

    if [ $? -ne 0 ]; then
	exit
    fi
    make -j ${PARALLEL_NP}
    if [ $? -ne 0 ]; then
	exit
    fi
    make install

    # Copy the library
    cp libmstk.a ${MSTK_PREFIX}/lib

}

################################################################################
#
# trilinos
#
################################################################################
function install_trilinos {
    TRILINOS_SRC=${PREFIX}/trilinos/trilinos-${TRILINOS_VERSION}-Source
    TRILINOS_BUILD=${PREFIX}/trilinos/trilinos-${TRILINOS_VERSION}-build
    TRILINOS_INSTALL=${PREFIX}/trilinos/trilinos-${TRILINOS_VERSION}-install

    rm -rf ${TRILINOS_BUILD} ${TRILINOS_INSTALL} ${TRILINOS_SRC}
    mkdir -p ${PREFIX}/trilinos
    tar ${TAR_FLAGS} ${SOURCE}/trilinos-${TRILINOS_VERSION}-Source.tar.gz -C ${PREFIX}/trilinos

    mkdir -p ${TRILINOS_BUILD}
    cd ${TRILINOS_BUILD}

    #export BOOST_ROOT=${BOOST_PREFIX}
#        -D Boost_INCLUDE_DIRS:FILEPATH=${BOOST_PREFIX}/include \
#        -D Boost_LIBRARY_DIRS:FILEPATH=${BOOST_PREFIX}/lib \

    cmake \
	-D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
	-D Trilinos_ENABLE_Fortran:BOOL=OFF \
	-D Trilions_ENABLE_TESTS:BOOL=OFF \
	-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
	-D Trilinos_ENABLE_Teuchos:BOOL=ON \
	-D Trilinos_ENABLE_Epetra:BOOL=ON \
	-D Trilinos_ENABLE_NOX:BOOL=ON \
	-D Trilinos_ENABLE_STK:BOOL=ON \
	-D Trilinos_ENABLE_Rythmos:BOOL=ON \
	-D Trilinos_ENABLE_Boost:BOOL=ON \
	-D BOOST_ROOT:FILEPATH=${BOOST_PREFIX} \
	-D NOX_ENABLE_EXAMPLES:BOOL=ON \
	-D Didasko_ENABLE_EXAMPLES:BOOL=OFF \
	-D TPL_ENABLE_MPI:BOOL=ON \
	-D MPI_BASE_DIR:PATH=${MPI_PREFIX} \
	-D TPL_ENABLE_ExodusII:BOOL=ON \
	-D TPL_ExodusII_LIBRARIES:FILEPATH=${PREFIX}/lib \
	-D TPL_ExodusII_INCLUDE_DIRS:FILEPATH=${PREFIX}/include \
	-D DART_TESTING_TIMEOUT:STRING=600 \
	-D CMAKE_BUILD_TYPE:STRING=DEBUG \
	-D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
	-D CMAKE_INSTALL_PREFIX:PATH=${TRILINOS_INSTALL} \
	${TRILINOS_SRC}

    if [ $? -ne 0 ]; then
	exit
    fi
    make -j${PARALLEL_NP}
    if [ $? -ne 0 ]; then
	exit
    fi
    ctest -W 100
    if [ $? -ne 0 ]; then
	exit
    fi
    make install
}

################################################################################
#
# ccse
#
################################################################################
function install_ccse {
    CCSE_DIR=${CCSE_PREFIX}/ccse/ccse-${CCSE_VERSION}
    CCSE_CONFIG=1
    CCSE_MAKE=1
    CCSE_TEST=1
    CCSE_INSTALL=1
    CCSE_VERBOSE=0
    CCSE_INSTALL_DIR=${PREFIX}/ccse/install

    rm -rf ${CCSE_DIR}
    mkdir -p ${CCSE_DIR}
    cd ${CCSE_PREFIX}/ccse
    tar zxf ${SOURCE}/ccse-${CCSE_VERSION}.tar.gz
    cd ${CCSE_DIR}

    # 1, 2, 3
    CCSE_SPACEDIM=${AMANZI_SPACEDIM}

    # FLOAT, DOUBLE
    CCSE_PRECISION=${AMANZI_PRECISION}

    if [ ${ENABLE_MPI} -eq 1 ]; then
	export MPI_EXEC=${MPI_PREFIX}/bin/mpiexec
	export MPI_EXEC_NUMPROCS_FLAG=-np
    fi

    mkdir -p ${CCSE_DIR}/build
    cd ${CCSE_DIR}/build

    if [ ${CCSE_VERBOSE} -eq 1 ]; then
	VFLAG="--debug-output"
    fi

    cmake \
	-D ENABLE_Config_Report:BOOL=ON \
	-D MPI_PREFIX:FILEPATH=${MPI_PREFIX} \
	-D ENABLE_TESTS:BOOL=ON \
	-D ENABLE_MPI:BOOL=${ENABLE_MPI} \
	-D ENABLE_OpenMP:BOOL=${ENABLE_OpenMP} \
	-D BL_SPACEDIM:INT=${CCSE_SPACEDIM} \
	-D BL_PRECISION:STRING="${CCSE_PRECISION}" \
        -D CMAKE_INSTALL_PREFIX:FILEPATH="${CCSE_INSTALL_DIR}" \
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

    make install
    if [ $? -ne 0 ]; then
	exit 1
    fi

    ctest --timeout 60 --output-on-failure
    if [ $? -ne 0 ]; then
	exit 1
    fi

}


################################################################################
#
# write a script for configuring and building amanzi
#
################################################################################
function generate_amanzi_build_script {
    cd ${SCRIPT_DIR}
    cat > amanzi-build.sh <<EOF
#!/bin/bash

AMANZI_WORK_DIR=${WORK_DIR}
AMANZI_CONFIG=0
AMANZI_DIR=
AMANZI_MAKE=0
AMANZI_INSTALL=0
AMANZI_CLOBBER=0
AMANZI_TEST=0
AMANZI_MAKE_NP=${PARALLEL_NP}
AMANZI_MAKE_VERBOSE=0
PLATFORM=${PLATFORM}
ENABLE_Structured=${ENABLE_Structured}
ENABLE_Unstructured=${ENABLE_Unstructured}

ENABLE_MPI=${ENABLE_MPI}
ENABLE_OpenMP=${ENABLE_OpenMP}

# 1, 2, 3
AMANZI_SPACEDIM=${AMANZI_SPACEDIM}

# FLOAT, DOUBLE
AMANZI_PRECISION=${AMANZI_PRECISION}

# ugly centos buildbot hack
CMAKE=`which cmake`
if [ -z \${CMAKE} ]; then
    if [ -n ${CMAKE_EXE} ]; then
	CMAKE=${CMAKE_EXE}
	echo "Using \${CMAKE}"
    else
	echo "Could not find cmake in path, set CMAKE_EXE in the config file."
	exit
    fi
fi
CTEST=`which ctest`
if [ -z \${CTEST} ]; then
    if [ -n ${CTEST_EXE} ]; then
	CTEST=${CTEST_EXE}
	echo "Using \${CTEST}"
    else
	echo "Could not find ctest in path, set CTEST_EXE in the config file."
	exit
    fi
fi

function determine_amanzi_dir() {
    AMANZI_DIR=
    case \$1 in
	/*) AMANZI_DIR=\$1;; # absolute path specified
	.) AMANZI_DIR=\$PWD;; # current directory specified
	*) AMANZI_DIR=\${AMANZI_WORK_DIR}/\$1;; # some other repo name specified
    esac
    #echo "Amanzi directory: \$AMANZI_DIR"
}

while getopts "abcd:imntv" flag
do
  case \$flag in
    a) AMANZI_CONFIG=1; AMANZI_MAKE=1; AMANZI_TEST=1; AMANZI_INSTALL=1;;
    b) AMANZI_CLOBBER=1;;
    c) AMANZI_CONFIG=1;;
    d) determine_amanzi_dir \${OPTARG};;
    i) AMANZI_INSTALL=1;;
    m) AMANZI_MAKE=1;;
    n) AMANZI_MAKE_NP= \${OPTARG};;
    t) AMANZI_TEST=1;;
    v) AMANZI_MAKE_VERBOSE=1;;
  esac
done

echo AMANZI_DIR=\$AMANZI_DIR
echo AMANZI_CONFIG=\$AMANZI_CONFIG
echo AMANZI_MAKE=\$AMANZI_MAKE
echo AMANZI_TEST=\$AMANZI_TEST

if [ \${PLATFORM} == Fedora -o \${PLATFORM} == RHEL ]; then
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
    # explicitly specify src/*.cmake so we do not accidently remove *.cmake files from tools!
    rm -rf \\
	CMakeCache.txt \\
	\${AMANZI_DIR}/src/*.cmake \${AMANZI_DIR}/src/*/*.cmake \${AMANZI_DIR}/src/*/*/*.cmake \\
	CMakeFiles */CMakeFiles */*/CMakeFiles \\
	Makefile */Makefile */*/Makefile \\
	Testing \\
	*~ */*~ */*/*~ \\
	*.rej */*.rej */*/*.rej
fi


if [ \$AMANZI_CONFIG -eq 1 ]; then
    mkdir -p \${AMANZI_DIR}/build
    mkdir -p \${AMANZI_DIR}/install
    cd \${AMANZI_DIR}/build

    if [ ${ENABLE_MPI} -eq 1 ]; then
        export CXX=${MPI_PREFIX}/bin/mpicxx
        export CC=${MPI_PREFIX}/bin/mpicc
        export FC=${MPI_PREFIX}/bin/mpif90
        export F90=${MPI_PREFIX}/bin/mpif90
        export F77=${MPI_PREFIX}/bin/mpif77
    fi

    \${CMAKE} \\
	-D CMAKE_INSTALL_PREFIX:FILEPATH=\${AMANZI_DIR}/install \\
	-D ENABLE_Config_Report:BOOL=ON \\
	-D ENABLE_MPI:BOOL=\${ENABLE_MPI} \\
	-D MPI_EXEC:FILEPATH=${MPI_PREFIX}/bin/mpiexec \\
	-D MPI_EXEC_NUMPROCS_FLAG:STRING=-np \\
	-D MPI_EXEC_GLOBAL_ARGS:STRING="${MPIEXEC_ARGS}" \\
	-D ENABLE_TESTS:BOOL=ON \\
	-D BOOST_ROOT:FILEPATH=${BOOST_PREFIX} \\
	-D UnitTest_DIR:FILEPATH=${UNITTEST_PREFIX} \\
	-D HDF5_DIR:FILEPATH=${HDF5_PREFIX} \\
	-D NetCDF_DIR:FILEPATH=${NETCDF_PREFIX} \\
	-D ExodusII_DIR:FILEPATH=${PREFIX} \\
	-D ENABLE_MESH_MOAB:BOOL=ON \\
	-D MOAB_DIR:FILEPATH=${PREFIX} \\
	-D ENABLE_MESH_MSTK:BOOL=${USE_MSTK} \\
	-D MSTK_DIR:FILEPATH=${MSTK_PREFIX} \\
	-D METIS_DIR:FILEPATH=${METIS_PREFIX} \\
	-D ENABLE_STK_Mesh:BOOL=ON \\
	-D Trilinos_DIR:FILEPATH=${PREFIX}/trilinos/trilinos-${TRILINOS_VERSION}-install \\
	-D ENABLE_OpenMP:BOOL=\${ENABLE_OpenMP} \\
	-D CCSE_DIR:FILEPATH=${CCSE_PREFIX}/ccse/install \\
	-D AMANZI_SPACEDIM:INT=\${AMANZI_SPACEDIM} \\
	-D AMANZI_PRECISION:STRING="\${AMANZI_PRECISION}" \\
	-D ENABLE_Structured:BOOL=\${ENABLE_Structured} \\
	-D ENABLE_Unstructured:BOOL=\${ENABLE_Unstructured} \\
	-D ASCEMIO_DIR:FILEPATH=${PREFIX} \\
        -D TESTS_REQUIRE_FULLPATH:BOOL=TRUE \\
	..

    if [ \$? -ne 0 ]; then
	exit 1
    fi
fi

if [ \${AMANZI_MAKE_VERBOSE} -eq 1 ]; then
    AMANZI_VFLAG="VERBOSE=1"
else
    AMANZI_VLFAG=
fi

if [ \$AMANZI_MAKE -eq 1 ]; then
    cd \${AMANZI_DIR}/build
    make -j \${AMANZI_MAKE_NP} \${AMANZI_VFLAG}
    if [ \$? -ne 0 ]; then
	exit 1
    fi
fi

if [ \$AMANZI_TEST -eq 1 ]; then
    cd \${AMANZI_DIR}/build
    \${CTEST} --timeout 150 --output-on-failure
    if [ \$? -ne 0 ]; then
	exit 1
    fi
fi

if [ \$AMANZI_INSTALL -eq 1 ]; then
    cd \${AMANZI_DIR}/build
    make install
    if [ \$? -ne 0 ]; then
        exit 1
    fi
fi

EOF
    chmod 700 amanzi-build.sh
}

################################################################################
#
# main
#
################################################################################

check_platform

CONFIG_FILE=${HOME}/.amanzi-tpls
DOWNLOAD_ARCHIVES=0
BUILD_OPENMPI=0
BUILD_BOOST=0
BUILD_HDF5=0
BUILD_ASCEMIO=0
BUILD_NETCDF=0
BUILD_EXODUS=0
BUILD_MOAB=0
BUILD_METIS=0
BUILD_MSTK=0
BUILD_TRILINOS=0
BUILD_ZLIB=0
BUILD_CURL=0
BUILD_UNITTEST=0
BUILD_CCSE=0
BUILD_CHECK=0

ENABLE_MPI=1
ENABLE_OpenMP=1
AMANZI_SPACEDIM=2
AMANZI_PRECISION=DOUBLE
ENABLE_Structured=0
ENABLE_Unstructured=1

while getopts "abcdefghikmnop:stuwz" flag
do
  case $flag in
    a) BUILD_OPENMPI=1; BUILD_BOOST=1; BUILD_CURL=1; BUILD_ZLIB=1; BUILD_HDF5=1; BUILD_ASCEMIO=1; BUILD_NETCDF=1; BUILD_EXODUS=1; BUILD_MOAB=1; BUILD_METIS=1; BUILD_MSTK=1; BUILD_TRILINOS=1; BUILD_UNITTEST=1; BUILD_CCSE=1;;
    b) BUILD_BOOST=1;;
    c) BUILD_CHECK=1;;
    d) DOWNLOAD_ARCHIVES=1;;
    e) BUILD_EXODUS=1;;
    f) BUILD_CCSE=1;;
    h) BUILD_HDF5=1;;
    i) BUILD_ASCEMIO=1;;
    k) BUILD_MSTK=1;;
    m) BUILD_MOAB=1;;
    n) BUILD_NETCDF=1;;
    o) BUILD_OPENMPI=1;;
    p) CONFIG_FILE=$OPTARG;;
    s) BUILD_METIS=1;;
    t) BUILD_TRILINOS=1;;
    u) BUILD_UNITTEST=1;;
    w) BUILD_CURL=1;;
    z) BUILD_ZLIB=1;;
  esac
done

if [ -e $CONFIG_FILE  ]; then
    . $CONFIG_FILE
else
    # the config file doesn't exist....
    echo "Could not find a valid tpl configuration file."
    exit
fi

echo "Using config file: $CONFIG_FILE"
echo "BUILD_OPENMPI=$BUILD_OPENMPI   MPI_PREFIX=$MPI_PREFIX"
echo "BUILD_BOOST=$BUILD_BOOST     BOOST_PREFIX=$BOOST_PREFIX"
echo "BUILD_CURL=$BUILD_CURL      CURL_PREFIX=$CURL_PREFIX"
echo "BUILD_ZLIB=$BUILD_ZLIB      ZLIB_PREFIX=$ZLIB_PREFIX"
echo "BUILD_HDF5=$BUILD_HDF5      HDF5_PREFIX=$HDF5_PREFIX"
echo "BUILD_ASCEMIO=$BUILD_ASCEMIO      ASCEMIO_PREFIX=$PREFIX"
echo "BUILD_NETCDF=$BUILD_NETCDF    NETCDF_PREFIX=$NETCDF_PREFIX"
echo "BUILD_EXODUS=$BUILD_EXODUS    EXODUS_PREFIX=$PREFIX"
echo "BUILD_MOAB=$BUILD_MOAB      MOAB_PREFIX=$PREFIX"
echo "BUILD_METIS=$BUILD_METIS     METIS_PREFIX=$PREFIX"
echo "BUILD_MSTK=$BUILD_MSTK      MSTK_PREFIX=$PREFIX"
echo "BUILD_UNITTEST=$BUILD_UNITTEST  UNITTEST_PREFIX=$UNITTEST_PREFIX"
echo "BUILD_TRILINOS=$BUILD_TRILINOS  TRILINOS_PREFIX=$PREFIX"
echo "BUILD_CCSE=$BUILD_CCSE      CCSE_PREFIX=$PREFIX"
echo "ENABLE_MPI=$ENABLE_MPI"
echo "ENABLE_OpenMP=$ENABLE_OpenMP"
echo "AMANZI_SPACEDIM=$AMANZI_SPACEDIM"
echo "AMANZI_PRECISION=$AMANZI_PRECISION"
echo "ENABLE_Structured=$ENABLE_Structured"
echo "ENABLE_Unstructured=$ENABLE_Unstructured"

if [ ${ENABLE_Structured} == 0 -a ${ENABLE_Unstructured} == 0 ]; then
  echo "Must enable Structured or Unstructured.  Exiting..."
  exit 1
fi

if [ $DOWNLOAD_ARCHIVES -eq 1 ]; then
    download_archives
    cd ${SCRIPT_DIR}
fi

if [ ${PLATFORM} == Fedora -o ${PLATFORM} == RHEL ]; then
    LD_RUN_PATH=$LD_RUN_PATH:${PREFIX}/lib
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PREFIX}/lib
    export LD_RUN_PATH LD_LIBRARY_PATH
fi

#
# build non-mpi dependent libraries and mpi using system compilers.
# The system compilers may not be the ones we want, so if the user
# specifies USE_XXX, then we use that.
#
if [ ! -z "${USE_CC}" ] && [ -n "${USE_CC}" ]; then
    # echo "USE_CC exists and is not empty"
    CC=${USE_CC}
else
    CC=gcc
fi
if [ ! -z "${USE_CXX}" ] && [ -n "${USE_CXX}" ]; then
    CXX=${USE_CXX}
else
    CXX=g++
fi
if [ ! -z "${USE_FC}" ] && [ -n "${USE_FC}" ]; then
    FC=${USE_FC}
else
    FC=gfortran
fi
if [ ! -z "${USE_F77}" ] && [ -n "${USE_F77}" ]; then
    F77=${USE_F77}
else
    F77=gfortran
fi
if [ ! -z "${USE_F90}" ] && [ -n "${USE_F90}" ]; then
    F90=${USE_F90}
else
    F90=gfortran
fi

export CC CXX FC F77 F90

if [ $BUILD_CURL -eq 1 ]; then
    install_curl
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_ZLIB -eq 1 ]; then
    install_zlib
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_UNITTEST -eq 1 ]; then
    install_unittest
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_BOOST -eq 1 ]; then
    install_boost
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_OPENMPI -eq 1 ]; then
    install_mpi
    cd ${SCRIPT_DIR}
fi

#
#  build mpi dependent libraries using mpi compilers
#
if [ $ENABLE_MPI -eq 1 ]; then
    export CXX=${MPI_PREFIX}/bin/mpicxx
    export CC=${MPI_PREFIX}/bin/mpicc
    export FC=${MPI_PREFIX}/bin/mpif90
    export F77=${MPI_PREFIX}/bin/mpif77
    export F90=${MPI_PREFIX}/bin/mpif90
fi

if [ $BUILD_HDF5 -eq 1 ]; then
    install_hdf5
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_ASCEMIO -eq 1 ]; then
    install_ascemio
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_NETCDF -eq 1 ]; then
    install_netcdf
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_EXODUS -eq 1 ]; then
    install_exodus
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_MOAB -eq 1 ]; then
#    download_moab
    install_moab
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_METIS -eq 1 ]; then
    install_metis
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_TRILINOS -eq 1 ]; then
    install_trilinos
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_MSTK -eq 1 ]; then
    install_mstk
    cd ${SCRIPT_DIR}
fi

if [ $BUILD_CCSE -eq 1 ]; then
    install_ccse
    cd ${SCRIPT_DIR}
fi

generate_amanzi_build_script

