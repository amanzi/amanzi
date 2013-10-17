source ~/.config_files/ats-modules

export CXX=`which mpiCC`
export CC=`which mpicc`
export BOOST_ROOT=/packages/boost/boost-1.46.1/openmpi-1.4.4_gcc-4.6.1

AMANZI_INSTALL_DIR=${AMANZI_DIR}/build-debug

ats_build=build-debug
rm -rf $ats_build
mkdir $ats_build
cd $ats_build

cmake \
    -D CMAKE_BUILD_TYPE:STRING=DEBUG \
    -D Amanzi_DIR=${AMANZI_INSTALL_DIR}/lib  \
    -D ATS_INSTALL_DIR=${ATS_DIR}/${ats_build} \
    -D MPI_DIR:FILEPATH=${MPI_HOME} \
    ..
