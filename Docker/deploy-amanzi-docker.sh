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
    echo "  --amanzi_tpls_ver   Which version of the Amanzi TPLs should we build? (Default: 0.98.9)"
    echo "  --output_style      Should we use the condensed or plain version of Docker output (Default: condensed)."
    echo "                      Set --output_style='plain' for expanded output"
    echo "  --multiarch         Build for both linux/amd64 and linux/arm64 instead of only local system architecture"
    echo "                      Assumes your already have Docker configured to build multiarchitecture images"
    echo "  --push              Push resulting image to Dockerhub (Requires caution!!)"
    echo "  --mpi_flavor        Which MPI implementation does the TPL image you want to build from use? (Default: mpich)"
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
    --amanzi_tpls_ver=*)
    amanzi_tpls_ver="${i#*=}"
    shift
    ;;
    --use_proxy=*)
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

AMANZI_GIT_LATEST_TAG_VER=`(cd $amanzi_src_dir; git tag -l amanzi-* | tail -n1 | sed -e 's/amanzi-//')`
AMANZI_GIT_GLOBAL_HASH=`(cd $amanzi_src_dir; git rev-parse --short HEAD)`
AMANZI_VER="${AMANZI_GIT_LATEST_TAG_VER}_${AMANZI_GIT_GLOBAL_HASH}"

echo ""
echo "AMANZI_SOURCE_DIR = $amanzi_src_dir"
echo " - latest tag       $AMANZI_GIT_LATEST_TAG_VER"
echo " - global hash      $AMANZI_GIT_GLOBAL_HASH"
echo " - version string   $AMANZI_VER"
echo ""

if $multiarch
then
    docker buildx build \
        --platform=linux/amd64,linux/arm64 \
        ${cache} \
        --build-arg amanzi_branch=${amanzi_branch} \
        --build-arg amanzi_tpls_ver=${amanzi_tpls_ver} \
        --build-arg mpi_flavor=${mpi_flavor} \
        ${use_proxy} \
        ${output} \
        ${push_arg} \
        -f ${amanzi_src_dir}/Docker/Dockerfile-Amanzi \
        -t metsi/amanzi:${AMANZI_VER} .
else
    docker build \
        ${cache} \
        --build-arg amanzi_branch=${amanzi_branch} \
        --build-arg amanzi_tpls_ver=${amanzi_tpls_ver} \
        --build-arg mpi_flavor=${mpi_flavor} \
        ${use_proxy} \
        ${output} \
        ${push_arg} \
        -f ${amanzi_src_dir}/Docker/Dockerfile-Amanzi \
        -t metsi/amanzi:${AMANZI_VER} .
fi

docker tag metsi/amanzi:${AMANZI_VER} metsi/amanzi:latest

