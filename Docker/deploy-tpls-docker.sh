#!/bin/bash

# Print out options and help statements
Help()
{
    echo "Usage: $(basename $0) [-h|--help] [additional options]"
    echo "Options:"
    echo "  -h, --help          Display this help message"
    echo "  --no_cache          Ignore docker layers in Docker build cache and build from scratch"
    echo "  --build_mpi         Build MPI implementation (either MPICH or OpenMPI) instead"
    echo "                      of using precompiled binaries in Ubuntu package repository (Default: True)"
    echo "  --mpi_distro        Which MPI implementation to use? (currently OpenMPI or MPICH) (Default: MPICH)"
    echo "  --mpi_version       Which version number of --mpi_distro to build (Default: 4.0.3)"
    echo "  --petsc_ver         Which version of PETSc to build as part of TPLs? (Default: 3.20)"
    echo "  --trilinos_ver      Which version of trilinos to build as part of TPLs? (Default: 15-1-6af5f44)"
    echo "  --amanzi_branch     Which Amanzi branch should be used when building container? (Default: master)"
    echo "  --amanzi_src_dir    Where does the Amanzi repo reside on the current system?"
    echo "                      (Default: /ascem/amanzi/repos/amanzi-master)"
    echo "  --amanzi_tpls_ver   Which version of the Amanzi TPLs should we build? (Default: 0.98.9)"
    echo "  --output_style      Should we use the condensed or plain version of Docker output (Default: condensed)"
    echo "  --multiarch         Build for both linux/amd64 and linux/arm64 instead of only local system architecture"
    echo "                      Assumes your already have Docker configured to build multiarchitecture images"
    echo "  --push              Push resulting image to Dockerhub (Requires caution!!)"
    exit 0
}

get_tpl_version()
#
#  Find version information and provide it for production installs
#
# rf notes 240528:
# are --petsc_ver and --trilinos_ver options still necessary?
# seems like these could be scraped from the TPL versions information,
# as in deploy-ats-docker.sh
{
   local SBFile=${amanzi_src_dir}/config/SuperBuild/TPLVersions.cmake
   tpl_version_major=`grep AMANZI_TPLS_VERSION_MAJOR ${SBFile} | tr -cd '[[:digit:]]'`
   tpl_version_minor=`grep AMANZI_TPLS_VERSION_MINOR ${SBFile} | tr -cd '[[:digit:]]'`
   tpl_version_patch=`grep AMANZI_TPLS_VERSION_PATCH ${SBFile} | tr -cd '[[:digit:]]'`
   echo "${tpl_version_major}.${tpl_version_minor}.${tpl_version_patch}"
}


# parse command line options, if given
for i in "$@"
do
case $i in
    -h|--help)
    Help
    ;;
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
    --amanzi_src_dir=*)
    amanzi_src_dir="${i#*=}"
    shift
    ;;
    --amanzi_tpls_ver=*)
    amanzi_tpls_ver="${i#*=}"
    shift
    ;;
    --use_proxy)
    use_proxy="--build-arg http_proxy=proxyout.lanl.gov:8080 --build-arg https_proxy=proxyout.lanl.gov:8080"
    shift
    ;;
    --output_style=*)
    output_style="${i#*=}"
    shift
    ;;
    --multiarch)
    multiarch=True
    shift
    ;;
    --push)
    push=True
    shift
    ;;
    --no_cache)
    cache="--no-cache"
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
amanzi_src_dir="${amanzi_src_dir:-/ascem/amanzi/repos/amanzi-master}"
amanzi_tpls_ver="${amanzi_tpls_ver:-`get_tpl_version`}"
use_proxy="${use_proxy:-}"
push="${push:-False}"
output_style="${output_style:-}"
multiarch="${multiarch:-False}"
cache="${cache:-}"

if [ "${output_style}" = "plain" ]; then
    output="--progress=plain"
else
    output=""
fi

if "${push}" ; then
    push_arg="--push"
else
    if $multiarch ; then
        push_arg="--load"
    else
        push_arg=""
    fi
fi

if $multiarch
then
    docker buildx build \
        --platform=linux/amd64,linux/arm64 \
        ${cache} \
        ${push_arg} \
        --build-arg petsc_ver=${petsc_ver} \
        --build-arg trilinos_ver=${trilinos_ver} \
        --build-arg amanzi_branch=${amanzi_branch} \
        --build-arg build_mpi=${build_mpi} \
        --build-arg mpi_flavor=${mpi_distro} \
        --build-arg mpi_version=${mpi_version} \
        ${use_proxy} \
        ${output} \
        -f ${amanzi_src_dir}/Docker/Dockerfile-TPLs \
        -t metsi/amanzi-tpls:${amanzi_tpls_ver}-${mpi_distro} .
else
    docker build \
        ${cache} \
        ${push_arg} \
        --build-arg petsc_ver=${petsc_ver} \
        --build-arg trilinos_ver=${trilinos_ver} \
        --build-arg amanzi_branch=${amanzi_branch} \
        --build-arg build_mpi=${build_mpi} \
        --build-arg mpi_flavor=${mpi_distro} \
        --build-arg mpi_version=${mpi_version} \
        ${use_proxy} \
        ${output} \
        -f ${amanzi_src_dir}/Docker/Dockerfile-TPLs \
        -t metsi/amanzi-tpls:${amanzi_tpls_ver}-${mpi_distro} .
fi
