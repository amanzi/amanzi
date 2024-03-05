#!/bin/bash

#
#  Set TPLs version for Docker image in Travis
#

show_usage()
{
cat << EOF
Usage: ${0##*/} [OPTIONS...]

Extract the Amanzi TPLs Collection version number,

    --amanzi-source-dir     - sets the amanzi repo directory (default=/home/metsi/amanzi)
    --help                  - Display this help and exit.

Or you can specify the Amanzi source directory in the AMANZI_SOURCE_DIR environment varaible.

EOF
}

#
# Default location based on Docker working-directory (relative path is ok).
#
AMANZI_SOURCE_DIR=../

parse_options()
{
    option_list=(help amanzi_source_dir:)

    params="$(getopt -o h --long "$(printf "%s," "${option_list[@]}")" --name "$0" -- "$@")"
    eval set -- "$params"

    while [[ $# -gt 0 ]] ; do
        case $1 in
            -h|--help)
                show_usage
                exit 0
                ;;
            --amanzi_source_dir)
                AMANZI_SOURCE_DIR=$2
                shift 2
                ;;
            --)
                shift
                break
                ;;
            *)
                echo -e "Error: $0 invalid option '$1'\nTry '$0 --help' for more information.\n" >&2
                exit 1
                ;;
        esac
        
    done
}


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

# Parse arguments
parse_options "$@"

# Get version
get_tpl_version

