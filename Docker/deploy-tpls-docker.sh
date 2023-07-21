#!/bin/bash

# MPI installed in the Docker image
# Options: openmpi, mpich
MPI_DISTRO=mpich

PETSC_VER=3.13
TRILINOS_VER=13-0-afc4e525

AMANZI_BRANCH=rfiorella/docker-multiarch
AMANZI_SOURCE_DIR=/amanzi/repos/amanzi-docker
AMANZI_TPLS_VER=0.98.6

# LANL_PROXY="--build-arg http_proxy=proxyout.lanl.gov:8080 --build-arg https_proxy=proxyout.lanl.gov:8080"

# --progress=plain

docker buildx build --no-cache --platform linux/arm64,linux/amd64 --build-arg petsc_ver=${PETSC_VER} --build-arg trilinos_ver=${TRILINOS_VER} --build-arg amanzi_branch=${AMANZI_BRANCH} -f ${AMANZI_SOURCE_DIR}/Docker/Dockerfile-TPLs -t metsi/amanzi-tpls:${AMANZI_TPLS_VER}-${MPI_DISTRO} --output type=registry .

# docker build --build-arg petsc_ver=${PETSC_VER} --build-arg trilinos_ver=${TRILINOS_VER} --build-arg amanzi_branch=${AMANZI_BRANCH} -f ${AMANZI_SOURCE_DIR}/Docker/Dockerfile-TPLs-base -t metsi/amanzi-tpls:${AMANZI_TPLS_VER}-${MPI_DISTRO} .
