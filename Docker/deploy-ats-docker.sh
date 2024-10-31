#!/bin/bash

# Print out options and help statements
Help()
{
    echo "Usage: $(basename $0) [-h|--help] [additional options]"
    echo "Options:"
    echo "  -h, --help          Display this help message"
    echo "  --no_cache          Ignore docker layers in Docker build cache and build from scratch"
    echo "  --amanzi_branch     Which Amanzi branch should be used when building container? (Default: master)"
    echo "  --amanzi_src_dir    Where does the Amanzi repo reside on the current system?"
    echo "                      (Default: /ascem/amanzi/repos/amanzi-master)"
    echo "  --ats_src_dir       Where does the ATS repo reside on the current system?"
    echo "                      (Default: \$amanzi_src_dir/src/physics/ats)"
    echo "  --amanzi_tpls_ver   Which version of the Amanzi TPLs should we build? (Default: 0.98.9)"
    echo "  --output_style      Should we use the condensed or plain version of Docker output (Default: condensed)."
    echo "                      Set --output_style='plain' for expanded output"
    echo "  --multiarch         Build for both linux/amd64 and linux/arm64 instead of only local system architecture"
    echo "                      Assumes your already have Docker configured to build multiarchitecture images"
    echo "  --mpi_flavor        Which MPI implementation does the TPL image you want to build from use? (Default: mpich)"
    echo "  --push              Push resulting image to Dockerhub (Requires caution!!)"
    echo "  --use_proxy         Passes https_proxy and http_proxy environment variables to docker build"
    exit 0
}

get_tpl_version()
#
#  Find version information and provide it for production installs
#
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
    --amanzi_branch=*)
    amanzi_branch="${i#*=}"
    shift
    ;;
    --amanzi_src_dir=*)
    amanzi_src_dir="${i#*=}"
    shift
    ;;
    --ats_src_dir=*)
    ats_src_dir="${i#*=}"
    shift
    ;;
    --amanzi_tpls_ver=*)
    amanzi_tpls_ver="${i#*=}"
    shift
    ;;
    --use_proxy)
    use_proxy="--build-arg http_proxy=${http_proxy} --build-arg https_proxy=${https_proxy}"
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
    --mpi_flavor=*)
    mpi_flavor="${i#*=}"
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
amanzi_branch="${amanzi_branch:-master}"
amanzi_src_dir="${amanzi_src_dir:-/ascem/amanzi/repos/amanzi-master}"
amanzi_tpls_ver="${amanzi_tpls_ver:-`get_tpl_version`}"
use_proxy="${use_proxy:-}"
output_style="${output_style:-}"
multiarch="${multiarch:-False}"
push="${push:-False}"
mpi_flavor="${mpi_flavor:-mpich}"
ats_src_dir="${ats_src_dir:-$amanzi_src_dir/src/physics/ats}"
cache="${cache:-}"

AMANZI_GIT_LATEST_TAG_VER=`(cd $amanzi_src_dir; git tag -l amanzi-* | tail -n1 | sed -e 's/amanzi-//')`
AMANZI_GIT_GLOBAL_HASH=`(cd $amanzi_src_dir; git rev-parse --short HEAD)`
AMANZI_VER="${AMANZI_GIT_LATEST_TAG_VER}_${AMANZI_GIT_GLOBAL_HASH}"

echo ""
echo "AMANZI_SRC_DIR    = $amanzi_src_dir"
echo " - latest tag       $AMANZI_GIT_LATEST_TAG_VER"
echo " - global hash      $AMANZI_GIT_GLOBAL_HASH"
echo " - version string   $AMANZI_VER"
echo ""

ATS_GIT_LATEST_TAG_VER=`(cd $ats_src_dir; git tag -l ats-* | tail -n1 | sed -e 's/ats-//')`
ATS_GIT_GLOBAL_HASH=`(cd $ats_src_dir; git rev-parse --short HEAD)`
ATS_VER="${ATS_GIT_LATEST_TAG_VER}_${ATS_GIT_GLOBAL_HASH}"

echo "ATS_SOURCE_DIR =    $ats_src_dir"
echo " - latest tag       $ATS_GIT_LATEST_TAG_VER"
echo " - global hash      $ATS_GIT_GLOBAL_HASH"
echo " - version string   $ATS_VER"
echo ""

if [ "${use_proxy}" ]; then
    LANL_PROXY="--build-arg http_proxy=proxyout.lanl.gov:8080 --build-arg https_proxy=proxyout.lanl.gov:8080"
else
    LANL_PROXY=""
fi
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
        ${output} \
        ${use_proxy} \
        ${push_arg} \
        --build-arg amanzi_branch=${amanzi_branch} \
        --build-arg amanzi_tpls_ver=${amanzi_tpls_ver} \
        --build-arg mpi_flavor=${mpi_flavor} \
        -f ${amanzi_src_dir}/Docker/Dockerfile-ATS-pyvista2 \
        -t metsi/ats:${ATS_VER} .
else
   docker build \
        ${cache} \
        ${output} \
        ${use_proxy} \
        ${push_arg} \
        --build-arg amanzi_branch=${amanzi_branch} \
        --build-arg amanzi_tpls_ver=${amanzi_tpls_ver} \
        --build-arg mpi_flavor=${mpi_flavor} \
        -f ${amanzi_src_dir}/Docker/Dockerfile-ATS-pyvista2 \
        -t metsi/ats:${ATS_VER} .
fi

docker tag metsi/ats:${ATS_VER} metsi/ats:latest

