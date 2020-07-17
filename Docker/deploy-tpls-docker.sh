#!/bin/bash

# MPI installed in the Docker image
# Options: openmpi, mpich
MPI_DISTRO=openmpi

PETSC_VER=3.11.3
TRILINOS_VER=12-18-55a7599733-Nov11

AMANZI_BRANCH=master
AMANZI_SOURCE_DIR=/ascem/amanzi/repos/amanzi-master
AMANZI_TPLS_VER=0.97.6

LANL_PROXY="--build-arg http_proxy=proxyout.lanl.gov:8080 --build-arg https_proxy=proxyout.lanl.gov:8080"

docker build --no-cache --build-arg petsc_ver=${PETSC_VER} --build-arg trilinos_ver=${TRILINOS_VER} --build-arg amanzi_branch=${AMANZI_BRANCH} -f ${AMANZI_SOURCE_DIR}/Docker/Dockerfile-TPLs -t metsi/amanzi-tpls:${AMANZI_TPLS_VER}-${MPI_DISTRO} .
