#!/bin/bash

# Print out options and help statements
Help()
{
    echo "Usage: $(basename $0) [-h|--help] [additional options]"
    echo "Options:"
    echo "  -h, --help          Display this help message"
    echo "  --amanzi_branch     Which Amanzi branch should be used when building container? (Default: master)"
    echo "  --amanzi_src_dir    Where does the Amanzi repo reside on the current system?"
    echo "                      (Default: /ascem/amanzi/repos/amanzi-master)"
    echo "  --amanzi_tpls_ver   Which version of the Amanzi TPLs should we build? (Default: 0.98.9)"
    echo "  --output_style      Should we use the condensed or plain version of Docker output (Default: condensed)."
    echo "                      Set --output_style='plain' for expanded output"
    echo "  --multiarch         Build for both linux/amd64 and linux/arm64 instead of only local system architecture"
    exit 0
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
    use_proxy="${i#*=}"
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
    *)
        # unknown option?
    ;;
esac
done

# set defaults, if not given on CLI
amanzi_branch="${amanzi_branch:-master}"
amanzi_src_dir="${amanzi_src_dir:-/ascem/amanzi/repos/amanzi-master}"
amanzi_tpls_ver="${amanzi_tpls_ver:-0.98.9}"
use_proxy="${use_proxy:-False}"
output_style="${output_style:-}"
multiarch="${multiarch:-False}"

if "${use_proxy}" ; then
    LANL_PROXY="--build-arg http_proxy=proxyout.lanl.gov:8080 --build-arg https_proxy=proxyout.lanl.gov:8080"
else
    LANL_PROXY=""
fi
if [ "${output_style}" = "plain" ]; then
    OUTPUT_STYLE="--progress=plain"
else
    OUTPUT_STYLE=""
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
    docker buildx build --platform=linux/amd64,linux/arm64 \
        --no-cache \
        --build-arg amanzi_branch=${amanzi_branch} \
        --build-arg amanzi_tpls_ver=${amanzi_tpls_ver} \
        ${LANL_PROXY} \
        ${OUTPUT_STYLE} \
        -f ${amanzi_src_dir}/Docker/Dockerfile-Amanzi \
        -t metsi/amanzi:${AMANZI_VER} .
else
    docker build --no-cache \
        --build-arg amanzi_branch=${amanzi_branch} \
        --build-arg amanzi_tpls_ver=${amanzi_tpls_ver} \
        ${LANL_PROXY} \
        ${OUTPUT_STYLE} \
        -f ${amanzi_src_dir}/Docker/Dockerfile-Amanzi \
        -t metsi/amanzi:${AMANZI_VER} .
fi

docker tag metsi/amanzi:${AMANZI_VER} metsi/amanzi:latest

