#!/bin/bash

# MPI installed in the Docker image
# Options: openmpi, mpich
MPI_DISTRO=mpich

PETSC_VER=3.13
TRILINOS_VER=13-0-afc4e525

AMANZI_BRANCH=master
AMANZI_SOURCE_DIR=/ascem/amanzi/repos/amanzi-master
AMANZI_TPLS_VER=0.97.12

LANL_PROXY="--build-arg http_proxy=proxyout.lanl.gov:8080 --build-arg https_proxy=proxyout.lanl.gov:8080"

docker build --no-cache --build-arg petsc_ver=${PETSC_VER} --build-arg trilinos_ver=${TRILINOS_VER} --build-arg amanzi_branch=${AMANZI_BRANCH} -f ${AMANZI_SOURCE_DIR}/Docker/Dockerfile-TPLs -t metsi/amanzi-tpls:${AMANZI_TPLS_VER}-${MPI_DISTRO} .
