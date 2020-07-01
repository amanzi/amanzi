#!/bin/bash

# PETSC_VER=3.11.3
# TRILINOS_VER=12-18-55a7599733-Nov11

get_tpl_version()
#
#  Find version information and provide it for production installs
#
{
   local SBFile=${AMANZI_SOURCE_DIR}/config/SuperBuild/TPLVersions.cmake
   tpl_version_major=`grep AMANZI_TPLS_VERSION_MAJOR ${SBFile} | tr -cd '[[:digit:]]'`
   tpl_version_minor=`grep AMANZI_TPLS_VERSION_MINOR ${SBFile} | tr -cd '[[:digit:]]'`
   tpl_version_patch=`grep AMANZI_TPLS_VERSION_PATCH ${SBFile} | tr -cd '[[:digit:]]'`
   echo "${tpl_version_major}.${tpl_version_minor}.${tpl_version_patch}"
}

AMANZI_BRANCH=master
AMANZI_SOURCE_DIR=/ascem/amanzi/repos/amanzi-master
AMANZI_TPLS_VER=`get_tpl_version`

ATS_SOURCE_DIR=$AMANZI_SOURCE_DIR/src/physics/ats


LANL_PROXY="--build-arg http_proxy=proxyout.lanl.gov:8080 --build-arg https_proxy=proxyout.lanl.gov:8080"

AMANZI_GIT_LATEST_TAG_VER=`(cd $AMANZI_SOURCE_DIR; git tag -l amanzi-* | tail -n1 | sed -e 's/amanzi-//')`
AMANZI_GIT_GLOBAL_HASH=`(cd $AMANZI_SOURCE_DIR; git rev-parse --short HEAD)`
AMANZI_VER="${AMANZI_GIT_LATEST_TAG_VER}_${AMANZI_GIT_GLOBAL_HASH}"

echo ""
echo "AMANZI_SOURCE_DIR = $AMANZI_SOURCE_DIR"
echo " - latest tag       $AMANZI_GIT_LATEST_TAG_VER"
echo " - global hash      $AMANZI_GIT_GLOBAL_HASH"
echo " - version string   $AMANZI_VER"
echo ""

ATS_GIT_LATEST_TAG_VER=`(cd $ATS_SOURCE_DIR; git tag -l ats-* | tail -n1 | sed -e 's/ats-//')`
ATS_GIT_GLOBAL_HASH=`(cd $ATS_SOURCE_DIR; git rev-parse --short HEAD)`
ATS_VER="${ATS_GIT_LATEST_TAG_VER}_${ATS_GIT_GLOBAL_HASH}"

echo "ATS_SOURCE_DIR =    $ATS_SOURCE_DIR"
echo " - latest tag       $ATS_GIT_LATEST_TAG_VER"
echo " - global hash      $ATS_GIT_GLOBAL_HASH"
echo " - version string   $ATS_VER"
echo ""

# MPI installed in the Docker image
# Options: openmpi, mpich
# MPI_FLAVOR=openmpi
MPI_FLAVOR=mpich

docker build --build-arg amanzi_branch=${AMANZI_BRANCH} --build-arg amanzi_tpls_ver=${AMANZI_TPLS_VER} --build-arg mpi_flavor=${MPI_FLAVOR} -f ${AMANZI_SOURCE_DIR}/Docker/Dockerfile-ATS -t metsi/ats:${ATS_VER} .

docker tag metsi/ats:${ATS_VER} metsi/ats:latest

