#!/bin/bash

# PETSC_VER=3.11.3
# TRILINOS_VER=12-18-55a7599733-Nov11

AMANZI_BRANCH=master
AMANZI_SOURCE_DIR=/ascem/amanzi/repos/amanzi-master
AMANZI_TPLS_VER=0.97.4

ATS_SOURCE_DIR=$AMANZI_SOURCE_DIR/src/physics/ats


LANL_PROXY="--build-arg http_proxy=proxyout.lanl.gov:8080 --build-arg https_proxy=proxyout.lanl.gov:8080"

AMANZI_GIT_LATEST_TAG_VER=`(cd $AMANZI_SOURCE_DIR; git tag -l amanzi-* | tail -n1 | sed -e 's/amanzi-//')`
AMANZI_GIT_GLOBAL_HASH=`(cd $AMANZI_SOURCE_DIR; git rev-parse --short HEAD)`
AMANZI_VER="${AMANZI_GIT_LATEST_TAG_VER}_${AMANZI_GIT_GLOBAL_HASH}"

echo $AMANZI_GIT_LATEST_TAG_VER
echo $AMANZI_GIT_GLOBAL_HASH
echo $AMANZI_VER

ATS_GIT_LATEST_TAG_VER=`(cd $ATS_SOURCE_DIR; git tag -l ats-* | tail -n1 | sed -e 's/ats-//')`
ATS_GIT_GLOBAL_HASH=`(cd $ATS_SOURCE_DIR; git rev-parse --short HEAD)`
ATS_VER="${ATS_GIT_LATEST_TAG_VER}_${ATS_GIT_GLOBAL_HASH}"

echo $ATS_GIT_LATEST_TAG_VER
echo $ATS_GIT_GLOBAL_HASH
echo $ATS_VER

docker build --build-arg amanzi_branch=${AMANZI_BRANCH} --build-arg amanzi_tpls_ver=${AMANZI_TPLS_VER} -f ${AMANZI_SOURCE_DIR}/Docker/Dockerfile-ATS -t metsi/ats:${ATS_VER} .

docker tag metsi/ats:${ATS_VER} metsi/ats:latest

