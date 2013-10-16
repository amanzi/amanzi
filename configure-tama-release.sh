source ~/.config_files/ats-modules

export CXX=`which mpiCC`
export CC=`which mpicc`
export BOOST_ROOT=/packages/boost/boost-1.46.1/openmpi-1.4.4_gcc-4.6.1

AMANZI_INSTALL_DIR=${AMANZI_DIR}/build-release

ats_build=build-release
rm -rf $ats_build
mkdir $ats_build
cd $ats_build

cmake \
    -D CMAKE_BUILD_TYPE:STRING=RELEASE \
    -D Amanzi_DIR=${AMANZI_INSTALL_DIR}/lib  \
    -D ATS_INSTALL_DIR=${ATS_DIR}/${ats_build} \
    -D ENABLE_DBC:BOOL=OFF \
    -D MPI_DIR:FILEPATH=${MPI_HOME} \
    ..
