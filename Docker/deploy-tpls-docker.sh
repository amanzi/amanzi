#!/bin/bash

# parse command line options, if given
for i in "$@"
do 
case $i in
    --build_mpi=*)
    build_mpi="${i#*=}"
    shift
    ;;
    --mpi_distro=*)
    mpi_distro="${i#*=}"
    shift
    ;;
    --mpi_version=*)
    mpi_version="${i#*=}"
    shift
    ;;
    --petsc_ver=*)
    petsc_ver="${i#*=}"
    shift
    ;;
    --trilinos_ver=*)
    trilinos_ver="${i#*=}"
    shift
    ;;
    --amanzi_branch=*)
    amanzi_branch="${i#*=}"
    shift
    ;;
    --amanzi_source_dir=*)
    amanzi_source_dir="${i#*=}"
    shift
    ;;
    --amanzi_tpls_ver=*)
    amanzi_tpls_ver="${i#*=}"
    shift
    ;;
    --use_proxy=*)
    use_proxy="${i#*=}"
    shift
    ;;
    --progress_style=*)
    progress_style="${i#*=}"
    shift
    ;;
    --multiarch=*)
    multiarch="${i#*=}"
    shift
    ;;
    *)
        # unknown option?
    ;;
esac
done

# set defaults, if not given on CLI
build_mpi="${build_mpi:-True}"
mpi_distro="${mpi_distro:-mpich}"
mpi_version="${mpi_version:-4.0.3}"
petsc_ver="${petsc_ver:-3.20}"
trilinos_ver="${trilnos_ver:-15-1-6af5f44}"
amanzi_branch="${amanzi_branch:-master}"
amanzi_source_dir="${amanzi_source_dir:-/ascem/amanzi/repos/amanzi-master}"
amanzi_tpls_ver="${amanzi_tpls_ver:-0.98.9}"
use_proxy="${use_proxy:-False}"
progress_style="${progress_style:-}"
multiarch="${multiarch:-False}"

if "${use_proxy}" ; then
    LANL_PROXY="--build-arg http_proxy=proxyout.lanl.gov:8080 --build-arg https_proxy=proxyout.lanl.gov:8080"
else
    LANL_PROXY=""
fi
if [ "${progress_style}" = "plain" ]; then
    PROGRESS_STYLE="--progress=plain"
else
    PROGRESS_STYLE=""
fi

if $multiarch
then
    docker buildx build --platform=linux/amd64,linux/arm64 \
             --no-cache --build-arg petsc_ver=${petsc_ver} \
             --build-arg trilinos_ver=${trilinos_ver} \
             --build-arg amanzi_branch=${amanzi_branch} \
             --build-arg build_mpi=${build_mpi} --build-arg mpi_flavor=${mpi_distro} \
             --build-arg mpi_version=${mpi_version} \
             $LANL_PROXY \
             $PROGRESS_STYLE \
             -f ${amanzi_source_dir}/Docker/Dockerfile-TPLs \
             -t metsi/amanzi-tpls:${amanzi_tpls_ver}-${mpi_distro} .
else
    docker build --no-cache --build-arg petsc_ver=${petsc_ver} \
             --build-arg trilinos_ver=${trilinos_ver} \
             --build-arg amanzi_branch=${amanzi_branch} \
             --build-arg build_mpi=${build_mpi} --build-arg mpi_flavor=${mpi_distro} \
             --build-arg mpi_version=${mpi_version} \
             $LANL_PROXY \
             $PROGRESS_STYLE \
             -f ${amanzi_source_dir}/Docker/Dockerfile-TPLs \
             -t metsi/amanzi-tpls:${amanzi_tpls_ver}-${mpi_distro} .
fi
