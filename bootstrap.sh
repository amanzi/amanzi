#!/bin/bash

# ############################################################################ #
#                                                                              #
#   Amanzi Boost Script                                                        #
#                                                                              #
#      Script that builds Amanzi from scratch                                  #
#                                                                              #
# ############################################################################ #

# ---------------------------------------------------------------------------- #
# Initialization
# ---------------------------------------------------------------------------- #

# Logic parameters
TRUE=1
FALSE=0

# System information
system_name=`uname`
system_arch=`uname -m`

# Test script controls
print_exit=${FALSE}

# Known compiler lists
known_c_compilers="gcc icc clang"
known_cxx_compilers="g++ icpc clang++"
known_fortran_compilers="gfortran ifort"

# Directory information
start_directory=$PWD
amanzi_source_dir=$(cd $(dirname "$0")/;pwd)
ats_source_dir="${amanzi_source_dir}/src/physics/ats"

# ASCEM Web address
ascem_protocol=https
ascem_site='github.com/amanzi'
ascem_tpl_site="${ascem_site}/amanzi-tpls"

# Default root build and install prefix
dflt_build_prefix=`pwd`
dflt_install_prefix=`pwd`

# Source and build directories
amanzi_build_dir="${dflt_build_prefix}/build/amanzi"
amanzi_install_prefix="${dflt_install_prefix}/install/amanzi"

# Git 
git_binary=`which git`

# CURL
curl_binary=`which curl`

# CMake
cmake_binary=`which cmake`
ctest_binary=`which ctest`
cmake_version=3.11.4
cmake_url=https://cmake.org/files/v3.11
cmake_archive_file=cmake-${cmake_version}.tar.gz

# Build configuration
parallel_jobs=2
build_type=Release
trilinos_build_type=Release
tpls_build_type=Release

# Compiler definitions
build_c_compiler=
build_cxx_compiler=
build_fort_compiler=

# Compiler flags
build_c_flags=
build_cxx_flags=
build_fort_flags=
build_link_flags=

# Accelerators
epetra=${TRUE}
kokkos=${FALSE}
kokkos_cuda=${FALSE}
kokkos_openmp=${FALSE}

# MPI installation
tools_mpi=openmpi
mpi_root_dir=
mpi_exec_args=

# TPL (Third Party Libraries)

# Spack
Spack=$FALSE
Spack_binary=`which spack`
build_Spack=$FALSE

# xSDK installation (optional)
xsdk=$FALSE #default is to not use
xsdk_root_dir=

# Tools build parameters
tools_build_dir=${dflt_build_prefix}/build/tools
tools_download_dir=${tools_build_dir}/Downloads
tools_install_prefix=${dflt_install_prefix}/install/tools

# Point to a configuration file
tpl_config_file=

# TPL build parameters
tpl_build_dir=${dflt_build_prefix}/build/tpls
tpl_download_dir=${tpl_build_dir}/Downloads
tpl_install_prefix=${dflt_install_prefix}/install/tpls

# build stages (TPLs, TPLs + Amanzi, TPLs + Amanzi + UserGuide)
build_tpls=$TRUE
build_amanzi=$TRUE
build_user_guide=$FALSE

# Color output display
no_color=$FALSE

# Amanzi build configuration
# -- components
structured=$TRUE
unstructured=$TRUE
geochemistry=$TRUE
# -- mesh frameworks
#    stk framewotk was deprecated and removed
mstk_mesh=$TRUE
moab_mesh=$FALSE
# -- tools
amanzi_branch=
ats_branch=
ccse_tools=$FALSE
spacedim=2
test_suite=$FALSE
reg_tests=$FALSE
build_stage_1=${FALSE}
build_stage_2=${FALSE}
dry_run=${FALSE}
# -- tpls
alquimia=${FALSE}
crunchtope=${FALSE}
hypre=${TRUE}
superlu=${TRUE}
netcdf4=${TRUE}
petsc=${TRUE}
pflotran=${FALSE}
amanzi_physics=${TRUE}
ats_physics=${FALSE}
shared=${FALSE}
silo=${FALSE}

# -- arch sets machine-specific variables
amanzi_arch=
prefer_static=${TRUE}
exec_static=${FALSE}
arch_tpl_opts=
arch_amanzi_opts=


# ---------------------------------------------------------------------------- #
# Functions: basic messages and exit functions
# ---------------------------------------------------------------------------- #

function exit_now()
{
  exit $1
}

function status_message()
{

  local GREEN='32m'
  if [ "${no_color}" -eq "${TRUE}" ]; then
    echo "[`date`] $1"
  else
    echo -e "[`date`]\033[$GREEN $1\033[m"
  fi  
}

function error_message()
{
  local RED='31m'
  if [ "${no_color}" -eq "${TRUE}" ]; then
    echo "Amanzi Bootstrap ERROR: $1" 1>&2
  else  
    echo -e "Amanzi Bootstrap ERROR:\033[$RED $1\033[m" 1>&2 
  fi
}

function warn_message()
{
  local PINK='35m'
  if [ "${no_color}" -eq "${TRUE}" ]; then
    echo "Amanzi Bootstrap Warning: $1" 1>&2
  else  
    echo -e "Amanzi Bootstrap Warning:\033[$PINK $1\033[m" 1>&2
  fi
}

function fatal_message()
{
  error_message $1
  exit_now 10
}


# ---------------------------------------------------------------------------- #
# Functions: utilities
# ---------------------------------------------------------------------------- #

function make_fullpath()
{
  if echo $1 | grep "^\/" > /dev/null 2> /dev/null; then
    echo $1
  else
    echo "$start_directory/$1"
  fi
}

function mkdir_now()
{
  echo "making $1"
  if [ ! -e $1 ]; then
    mkdir -p $1
  fi
}

function download_file()
{
  url=$1
  file=$2
  output=$3

  curl_opts=
  if [ ! -z "${output}" ]; then
    if [ -d "${output}" ]; then
      curl_opts="--output $output/$file"
    else
      curl_opts="--output $output"
    fi
  else
    curl_opts="--remote-name"
  fi
  cmd="${curl_binary} ${curl_opts} $url/$file"
  echo $cmd
  ${curl_binary} $curl_opts $url/$file
  status=$?

  if [ "${status}" -ne "0" ]; then
     error_message "${curl_binary} returned $status"
     error_message "Failed to download $file from $url"
     exit_now $status
  fi
}


# ---------------------------------------------------------------------------- #
# Functions: command line functions
# ---------------------------------------------------------------------------- #

function set_feature()
{
  feature=$1
  action=$2

  if echo $action | grep "disable" > /dev/null 2>/dev/null; then
    eval "${feature}=$FALSE"
  fi
  if echo $action | grep "enable" > /dev/null 2>/dev/null; then
    eval "${feature}=$TRUE"
  fi
}

function parse_feature()
{
  feature_opt=$1

  if echo ${feature_opt} | grep "^--disable-" > /dev/null 2> /dev/null; then
    echo ${feature_opt} | sed "s/^--disable-//"
  fi
  if echo ${feature_opt} | grep "^--enable-" > /dev/null 2> /dev/null; then
    echo ${feature_opt} | sed "s/^--enable-//"
  fi
}

function parse_option_with_equal()
{
  a=$1
  opt_name=$2

  if echo $a | grep "^--${opt_name}=" > /dev/null 2> /dev/null; then
    echo $a | sed "s/^--${opt_name}=//"
  fi
}

function print_help()
{
echo '
Usage: '"$0"' [options]
Options: [defaults in brackets after descriptions]
Configuration:

  --help                  print this message
  
  --parallel=n            build in parallel, where n is
                          number of maximum make jobs ['"${parallel_jobs}"']

  --no-color              deactivate color-coded status output ['"${no_color}"']                        

  --tpl-config-file=FILE  define a CMake TPL configuration file. If this
                          option is selected, '"$0"' will NOT build the TPLs.

  --opt                   build optimized TPLs and Amanzi binaries. This the default
                          configuration.

  --relwithdebinfo        build optimized TPLs and Amanzi binaries with debug info
                          (for profiling)

  --debug                 build debug TPLs and Amanzi binaries.

  --branch=BRANCH         build TPLs and Amanzi found in BRANCH ['"${amanzi_branch}"']

  --branch_ats=BRANCH     build ATS found in BRANCH ['"${ats_branch}"']
  
  --spacedim=DIM          dimension of structured build (DIM=2 or 3) ['"${spacedim}"']

  --dry_run               show the configuration commands (but do not execute them) ['"${dry_run}"']

  --arch                  define the architecture for the build, only Summit, NERSC available


Build features:
Each feature listed here can be enabled/disabled with --[enable|disable]-[feature]
Value in brackets indicates default setting.

  structured              build structured AMR mesh capability ['"${structured}"']
  ccse_tools              build structured AMR tools for post processing and tecplot ['"${ccse_tools}"']

  unstructured            build unstructured mesh capability ['"${unstructured}"']
  mstk_mesh               build the MSTK Mesh Toolkit ['"${mstk_mesh}"']
  moab_mesh               build the MOAB Mesh Toolkit ['"${moab_mesh}"']

  hypre                   build the HYPRE solver APIs ['"${hypre}"']
  superlu                 build the SuperLU solver ['"${superlu}"']
  petsc                   build the PETSc solver APIs ['"${petsc}"']

  geochemistry            build all geochemistry packages (pflotran, crunchtope, alquimia) ['"${geochemistry}"']
  pflotran                build the PFlotran geochemistry backend ['"${pflotran}"']
  crunchtope              build the CrunchTope geochemistry backend ['"${crunchtope}"']
  alquimia                build the Alquimia geochemistry solver APIs ['"${alquimia}"']

  amanzi_physics          build Amanzi native physics package ['"${amanzi_physics}"']
  ats_physics             build ATS physics package (currently mutually exclusive) ['"${ats_physics}"']

  test_suite              run Amanzi Test Suite before installing ['"${test_suite}"']
  reg_tests               build regression tests into Amanzi Test Suite ['"${reg_tests}"']
  shared                  build Amanzi and tpls using shared libraries ['"${shared}"']
  Spack                   build TPLs using the Spack package manager when appropriate ['"${Spack}"']
  xsdk                    build TPLs available in xSDK first, then supplement with additional 
                          individual TPL builds ['"${xsdk}"']

  epetra                  build the Epetra stack of TPLs ['"${epetra}"']
  kokkos                  build the Tpetra/Kokkos stack of TPLs ['"${kokkos}"']
  kokkos_cuda             build Kokkos with Cuda node support (currently only Kokkos) ['"${kokkos_cuda}"']
  kokkos_openmp           build Kokkos with OpenMP node support (currently only Kokkos) ['"${kokkos_openmp}"']

  build_amanzi            build TPLs and Amanzi ['"${build_amanzi}"']
  build_user_guide        build TPLs, Amanzi, and UserGuide ['"${build_user_guide}"']


Tool definitions:

  --with-c-compiler=FILE     FILE is the C compiler
  --with-cxx-compiler=FILE   FILE is the C++ compiler
  --with-fort-compiler=FILE  FILE is the Fortran compiler

  --with-c-flags=STRING      STRING is additional C compiler flags
  --with-cxx-flags=STRING    STRING is additional C++ compiler flags
  --with-fort-flags=STRING   STRING is additional Fortran compiler flags

  --with-link-flags=STRING   STRING is additional linker flags

  --with-cmake[=FILE]        FILE is the CMake binary ['"${cmake_binary}"'] without FILE builds CMake
  --with-ctest=FILE          FILE is the CTest binary ['"${ctest_binary}"'], ignored if --with-cmake is set
  --with-git=FILE            FILE is the git binary ['"${git_binary}"']
  --with-curl=FILE           FILE is the CURL binary ['"${curl_binary}"']

  --with-mpi=DIR             use MPI installed in DIR ['"${mpi_root_dir}"'].
                             Will search for MPI compiler wrappers under this directory. It this
                             option is missing, OpenMPI will be built using the provided compiler.

  --with-xsdk=DIR            use libraries already available in xSDK installation in lieu of
                             downloading and installing them individually. ['"${xsdk_root_dir}"']


Directories and file names: 

  --prefix=PREFIX                install ALL files in tree rooted at PREFIX
                                 ['"${dflt_install_prefix}"']

  --amanzi-install-prefix=DIR    install Amanzi in tree rooted at DIR
                                 ['"${amanzi_install_prefix}"']

  --amanzi-build-dir=DIR         build Amanzi in DIR
                                 ['"${amanzi_build_dir}"']

  --tpl-install-prefix=DIR       install Amanzi TPLs in tree rooted at DIR
                                 ['"${tpl_install_prefix}"']
  
  --tpl-build-dir=DIR            build Amanzi TPLs in DIR
                                 ['"${tpl_build_dir}"']

  --tpl-download-dir=DIR         direct downloads of Amanzi TPLs to DIR
                                 ['"${tpl_download_dir}"']

  --tools-install-prefix=DIR     install Amanzi tools in tree rooted at DIR
                                 ['"${tools_install_prefix}"']. The list of tools includes
                                 currently OpenMPI and MPICH.
  
  --tools-build-dir=DIR          build Amanzi tools in DIR
                                 ['"${tools_build_dir}"']

  --tools-download-dir=DIR       direct downloads of Amanzi tools to DIR
                                 ['"${tools_download_dir}"']

  --tools-mpi=NAME               implementation of the Message Passing Interface (MPI)
                                 standard. NAME is either openmpi (default) or mpich.


Example with pre-existing MPI installation that builds dynamic libraries
and executables for external packages (TPLs) and then for Amanzi:

  ./bootstrap.sh --tpl-install-prefix=$HOME/TPLs-0.96.3-clang-9.0.0-ompi-3.1.4
                 --with-mpi=/opt/local 
                 --parallel=8 
                 --enable-shared
                 --enable-alquimia --enable-pflotran --enable-crunchtope 
                 --enable-petsc

Example that builds first OpenMPI, then TPLs, and finally Amanzi. It uses 
OSX C and C++ compilers and Fortran compiler from MacPorts:

  ./bootstrap.sh --tpl-install-prefix=$HOME/TPLs-0.96.3-clang-9.0.0-ompi-3.1.4
                 --with-c-compiler=/usr/bin/clang
                 --with-cxx-compiler=/usr/bin/clang++
                 --with-fort-compiler=/opt/local/bin/gfortran-mp-6
                 --parallel=8 
                 --enable-alquimia --enable-pflotran --enable-crunchtope 
                 --enable-petsc
                 --tools-mpi=openmpi
'
}

function print_variable_values
{
echo '

List of computed and ADJUSTED parameters
========================================
Compilers and Flags:
    build_c_compiler    = '"${build_c_compiler}"'
    build_cxx_compiler  = '"${build_cxx_compiler}"'
    build_fort_compiler = '"${build_fort_compiler}"'
    build_c_flags       = '"${build_c_flags}"'
    build_cxx_flags     = '"${build_cxx_flags}"'
    build_fort_flags    = '"${build_fort_flags}"'
    build_link_flags    = '"${build_link_flags}"'
    mpi_root_dir        = '"${mpi_root_dir}"'

Build configuration:
    build_type          = '"${build_type}"'
    build_stage_1       = '"${build_stage_1}"'
    build_stage_2       = '"${build_stage_2}"'
    parallel            = '"${parallel_jobs}"'
    shared              = '"${shared}"'
    static_libs_prefer  = '"${prefer_static}"'
    static_executables  = '"${exec_static}"'
    trilinos_build_type = '"${trilinos_build_type}"'
    tpls_build_type     = '"${tpls_build_type}"'
    tpl_config_file     = '"${tpl_config_file}"'
    amanzi_arch         ='"${amanzi_arch}"'

Amanzi Components:   
    structured     = '"${structured}"'
    spacedim       = '"${spacedim}"'
    unstructured   = '"${unstructured}"'
    amanzi_physics = '"${amanzi_physics}"'
    ats_physics    = '"${ats_physics}"'

Amanzi TPLs:
    alquimia     = '"${alquimia}"'
    crunchtope   = '"${crunchtope}"'
    geochemistry = '"${geochemistry}"'
    mstk_mesh    = '"${mstk_mesh}"'
    moab_mesh    = '"${moab_mesh}"'
    netcdf4      = '"${netcdf4}"'
    hypre        = '"${hypre}"'
    superlu      = '"${superlu}"'
    petsc        = '"${petsc}"'
    epetra       = '"${epetra}"'
    kokkos       = '"${kokkos}"' (Cuda='"${kokkos_cuda}"', OpenMP='"${kokkos_openmp}"')
    pflotran     = '"${pflotran}"'
    silo         = '"${silo}"'
    Spack        = '"${Spack}"'
    xsdk         = '"${xsdk}"'

Tools and Tests:
    ccse_tools   = '"${ccse_tools}"'
    cmake_binary = '"${cmake_binary}"'
    ctest_binary = '"${ctest_binary}"'
    curl_binary  = '"${curl_binary}"'
    git_binary   = '"${git_binary}"'
    Spack_binary = '"${Spack_binary}"'
    reg_tests    = '"${reg_tests}"'
    test_suite   = '"${test_suite}"'
    tools_mpi    = '"${tools_mpi}"'

Directories:
    prefix                = '"${prefix}"'
    amanzi_install_prefix = '"${amanzi_install_prefix}"'
    amanzi_build_dir      = '"${amanzi_build_dir}"'
    tpl_install_prefix    = '"${tpl_install_prefix}"'
    tpl_build_dir         = '"${tpl_build_dir}"'
    tpl_download_dir      = '"${tpl_download_dir}"'
    tools_install_prefix  = '"${tools_install_prefix}"'
    tools_build_dir       = '"${tools_build_dir}"'
    tools_download_dir    = '"${tpl_download_dir}"'
    
'    
}  

function parse_argv()
{
echo '
List of INPUT parameters
========================'
  argv=( "$@" )
  last=$(( ${#argv[@]} - 1 ))
  i=0
  while [ $i -le ${last} ]
  do
    opt=${argv[$i]}
    echo "opt:  $opt"
    case ${opt} in

      -h|--h|--help)
                print_help
                exit_now 0
                ;;

      --prefix=*)
                 tmp=`parse_option_with_equal "${opt}" 'prefix'`
                 prefix=`make_fullpath ${tmp}`
                 ;;

      --arch=*)
                 amanzi_arch=`parse_option_with_equal "${opt}" 'arch'`
		 ;;

      --parallel=[0-9]*)
                 parallel_jobs=`parse_option_with_equal "${opt}" 'parallel'`
                 ;;
      
      --opt)
                 build_type=Release
                 ;;
      --opt_tpls)
                 tpls_build_type=Release
                 ;;
      --opt_trilinos)
                 trilinos_build_type=Release
                 ;;

      --relwithdebinfo)
                 build_type=RelWithDebInfo
                 ;;
      --relwithdebinfo_tpls)
                 tpls_build_type=RelWithDebInfo
                 ;;
      --relwithdebinfo_trilinos)
                 trilinos_build_type=RelWithDebInfo
                 ;;
      --debug)
                 build_type=Debug
                 ;;
      --debug_tpls)
                 tpls_build_type=Debug
                 ;;
      --debug_trilinos)
                 trilinos_build_type=Debug
                 ;;
      --dry_run)
                 dry_run=${TRUE}
                 ;;

      --disable-*)
                 feature=`parse_feature "${opt}"`
                 set_feature ${feature} 'disable'
                 ;;

      --enable-*)
                 feature=`parse_feature "${opt}"`
                 set_feature ${feature} 'enable'
                 ;;

      --no-color)
                 no_color=${TRUE}
                 ;;

      --branch=*)
                 amanzi_branch=`parse_option_with_equal "${opt}" 'branch'`
                 ;;

      --branch_ats=*)
                 ats_branch=`parse_option_with_equal "${opt}" 'branch_ats'`
                 ;;

      --spacedim=*)
                 spacedim=`parse_option_with_equal "${opt}" 'spacedim'`
                 ;;

      --with-c-compiler=*)
                 tmp=`parse_option_with_equal "${opt}" 'with-c-compiler'`
                 build_c_compiler=`make_fullpath $tmp`
                 ;;

      --with-c-flags=*)
                 build_c_flags=`parse_option_with_equal "${opt}" 'with-c-flags'`
                 ;;

      --with-cxx-compiler=*)
                 tmp=`parse_option_with_equal "${opt}" 'with-cxx-compiler'`
                 build_cxx_compiler=`make_fullpath $tmp`
                 ;;

      --with-cxx-flags=*)
                 build_cxx_flags=`parse_option_with_equal "${opt}" 'with-cxx-flags'`
                 ;;

      --with-fort-compiler=*)
                 tmp=`parse_option_with_equal "${opt}" 'with-fort-compiler'`
                 build_fort_compiler=`make_fullpath $tmp`
                 ;;

      --with-fort-flags=*)
                 build_fort_flags=`parse_option_with_equal "${opt}" 'with-fort-flags'`
                 ;;

      --with-link-flags=*)
                 build_link_flags=`parse_option_with_equal "${opt}" 'with-link-flags'`
                 ;;

      --with-cmake=*)
                 tmp=`parse_option_with_equal "${opt}" 'with-cmake'`
                 cmake_binary=`make_fullpath $tmp`
                 ;;

      --with-ctest=*)
                 tmp=`parse_option_with_equal "${opt}" 'with-ctest'`
                 ctest_binary=`make_fullpath $tmp`
                 ;;

      --with-cmake)
                 cmake_binary=
                 ;;

      --with-git=*)
                 tmp=`parse_option_with_equal "${opt}" 'with-git'`
                 git_binary=`make_fullpath $tmp`
                 ;;

      --with-curl=*)
                 tmp=`parse_option_with_equal "${opt}" 'with-curl'`
                 curl_binary=`make_fullpath $tmp`
                 ;;

      --with-mpi=*)
                 tmp=`parse_option_with_equal "${opt}" 'with-mpi'`
                 mpi_root_dir=$tmp
                 ;;

      --with-xsdk=*)
                 tmp=`parse_option_with_equal "${opt}" 'with-xsdk'`
                 xsdk_root_dir=`make_fullpath $tmp`
                 XSDK=TRUE
                 ;;

      --amanzi-build-dir=*)
                 tmp=`parse_option_with_equal "${opt}" 'amanzi-build-dir'`
                 amanzi_build_dir=`make_fullpath $tmp`
                  ;;

      --amanzi-install-prefix=*)
                 tmp=`parse_option_with_equal "${opt}" 'amanzi-install-prefix'`
                 amanzi_install_prefix=`make_fullpath $tmp`
                 ;;

      --tpl-install-prefix=*)
                 tmp=`parse_option_with_equal "${opt}" 'tpl-install-prefix'`
                 tpl_install_prefix=`make_fullpath $tmp`
                 ;;

      --tpl-build-dir=*)
                 tmp=`parse_option_with_equal "${opt}" 'tpl-build-dir'`
                 tpl_build_dir=`make_fullpath $tmp`
                 ;;

      --tpl-download-dir=*)
                 tmp=`parse_option_with_equal "${opt}" 'tpl-download-dir'`
                 tpl_download_dir=`make_fullpath $tmp`
                 ;;
      
      --tpl-config-file=*)
                 tmp=`parse_option_with_equal "${opt}" 'tpl-config-file'`
                 tpl_config_file=`make_fullpath $tmp`
                 ;;

      --build_user_guide)
                 build_user_guide=${TRUE}
                 ;;

      --tools-install-prefix=*)
                 tmp=`parse_option_with_equal "${opt}" 'tools-install-prefix'`
                 tools_install_prefix=`make_fullpath $tmp`
                 ;;

      --tools-build-dir=*)
                 tmp=`parse_option_with_equal "${opt}" 'tools-build-dir'`
                 tools_build_dir=`make_fullpath $tmp`
                 ;;

      --tools-download-dir=*)
                 tmp=`parse_option_with_equal "${opt}" 'tools-download-dir'`
                 tools_download_dir=`make_fullpath $tmp`
                 ;;
      
      --tools-mpi=*)
                 tools_mpi=`parse_option_with_equal "${opt}" 'tools-mpi'`
                 ;;

      --print)
                 print_exit=${TRUE}
                 ;;

       *)
                 error_message "'${opt}' is an unknown option or an option missing a value"
                 exit_now 20
                 ;;
      esac

      i=$[$i+1]
  done

  # enforce implicit rules
  if [ "${build_user_guide}" -eq "${TRUE}" ]; then
    build_amanzi=${TRUE}
  fi

  if [ "${epetra}" -eq "${FALSE}" ]; then
      warning_message "Disabling Epetra disables all of geochemistry and hypre"
      geochemistry=${FALSE}
      hypre=${FALSE}
      superlu=${FALSE}
      alquimia=${FALSE}
      pflotran=${FALSE}
      crunchtope=${FALSE}
      petsc=${FALSE}
  fi
  
  if [ "${geochemistry}" -eq "${TRUE}" ]; then
    alquimia=${geochemistry}
    pflotran=${geochemistry}
    crunchtope=${geochemistry}
  fi

  # check compatibility
  if [ "${geochemistry}" -eq "${FALSE}" ]; then
    if [ "${pflotran}" -eq "${TRUE}" ]; then
      error_message "Option 'geochemisty' is incomatible with option 'pflotran'"
      exit_now 30
    fi
    if [ "${alquimia}" -eq "${TRUE}" ]; then
      error_message "Option 'geochemisty' is incomatible with option 'alquimia'"
      exit_now 30
    fi
    if [ "${crunchtope}" -eq "${TRUE}" ]; then
      error_message "Option 'geochemisty' is incomatible with option 'crunchtope'"
      exit_now 30
    fi
  fi

  # superlu required by petsc and hypre
  if [ "${petsc}" -eq "${TRUE}" ]; then
    if [ "${superlu}" -eq "${FALSE}" ]; then
      warning_message "Option 'petsc' requires 'superlu', turning it on."
      superlu=${TRUE}
    fi
  fi
  if [ "${hypre}" -eq "${TRUE}" ]; then
    if [ "${superlu}" -eq "${FALSE}" ]; then
      warning_message "Option 'hypre' requires 'superlu', turning it on."
      superlu=${TRUE}
    fi
  fi
  
  # deprecated options
  if [ "${tpls_only}" ]; then
    error_message "Option 'tpls_only' has been deprecated"
    exit_now 30
  fi
}


# ---------------------------------------------------------------------------- #
# CURL functions
# ---------------------------------------------------------------------------- #

curl_version()
{
  $curl_binary --version
}

curl_protocols()
{
  vers_string=`$curl_binary --version | grep Protocols: | sed 's/Protocols://'`

  echo ${vers_string}
}

curl_test_protocol()
{
  local protocols tmp type=$1
  protocols=`curl_protocols`
  tmp=`echo $protocols | sed "s/.*${type}.*/${type}-TRUE/g"`
  if [ "${tmp}" = "${type}-TRUE" ]; then
    echo ${TRUE}
  else
    echo ${FALSE}
  fi
}  


# ---------------------------------------------------------------------------- #
# CMake functions
# ---------------------------------------------------------------------------- #

function cmake_version
{
  ver_string=`${cmake_binary} --version | sed 's/cmake version[ ]*//'`
  echo ${ver_string}
}

function download_cmake
{
   target_dir=$1
   target=
   pwd_save=`pwd`
   if [ ! -z "${target_dir}" ]; then
     mkdir_now ${target_dir}
     cd ${target_dir}
     target=${target_dir}/${cmake_archive_file}
   fi

   download_file ${cmake_url} ${cmake_archive_file} ${target}

   cd ${pwd_save}

   status_message "CMake download complete"
}

function unpack_cmake
{
  archive=$1
  final_source_dir=$2

  pwd_save=`pwd`

  # Create the final source destination
  mkdir_now ${final_source_dir}
  cd ${final_source_dir}

  # Untar the archive file removing the leading path component  
  tar -z -x -f ${archive} --strip-components 1

  if [ $? -ne 0 ]; then
    error_message "Failed to untar CMake archive file ${archive}"
    exit_now 30
  fi

  cd ${pwd_save}
  status_message "Unpack CMake complete"
}

function configure_cmake
{
  source_dir=$1
  build_dir=$2
  install_dir=$3

  pwd_save=`pwd`

  # Create and cd to the build directory
  mkdir_now ${build_dir}
  cd ${build_dir}

  # Build the bootstrap args
  boot_args="--prefix=${install_dir}"
  echo "boot_args=$boot_args"

  curl_supports_https=`curl_test_protocol https`
  if [ "${curl_supports_https}" ]; then
    boot_args="${boot_args} --system-curl"
  fi

  boot_args="${boot_args} --parallel=${parallel_jobs}"

  # Run the bootstrap script
  bootstrap_script=${source_dir}/bootstrap
  if [ ! -e ${bootstrap_script} ]; then
    error_message "CMake bootstrap script ${bootstrap_script} does not exist"
    exit_now 30
  fi

  ${bootstrap_script} ${boot_args}

  if [ $? -ne 0 ]; then
    error_message "CMake bootstrap script failed"
    exit_now 30
  fi

  cd ${pwd_save}

  status_message "Configure CMake complete"
}

function build_cmake
{
  root_build_dir=$1
  root_install_dir=$2
  root_download_dir=$3

  cmake_root_build_dir=${root_build_dir}/cmake
  cmake_source_dir=${cmake_root_build_dir}/cmake-${cmake_version}-source
  cmake_build_dir=${cmake_root_build_dir}/cmake-${cmake_version}-build

  # Download the file
  download_file ${cmake_url} ${cmake_archive_file} ${root_download_dir}

  # Define the full path name of the archive file
  if [ ! -z "${root_download_dir}" ]; then
    archive=${root_download_dir}/${cmake_archive_file}
  else
    archive=${cmake_archive_file}
  fi

  # Unpack the archive file
  unpack_cmake ${archive} ${cmake_source_dir}

  # Configure (bootstrap) CMake
  configure_cmake ${cmake_source_dir} ${cmake_build_dir} ${root_install_dir}

  # Build CMake 
  cd ${cmake_build_dir}
  make -j ${parallel_jobs}
  if [ $? -ne 0 ]; then
    error_message "Failed to build CMake"
    exit_now 30
  fi

  # Install CMake 
  make install
  if [ $? -ne 0 ]; then
    error_message "Failed to install CMake"
    exit_now 30
  fi

  cmake_binary=${root_install_dir}/bin/cmake
  ctest_binary=${root_install_dir}/bin/ctest

  status_message "CMake build complete"
}


# ---------------------------------------------------------------------------- #
# Git functions
# ---------------------------------------------------------------------------- #

function ascem_git_clone
{
  repo=$1
  ${git_binary} clone ${ascem_protocol}://${ascem_site}/$repo
  if [ $? -ne 0 ]; then
    error_message "Failed to clone ${repo} from ${ascem_site}"
    exit_now 30
  fi
}

function git_change_branch()
{
  branch=$1
  save_dir=`pwd`
  cd ${amanzi_source_dir}
  status_message "In ${amanzi_source_dir} checking out ${branch}"
  ${git_binary} checkout ${branch}
  if [ $? -ne 0 ]; then
    error_message "Failed to update ${amanzi_source_dir} to branch ${branch}"
    exit_now 30
  fi
  cd ${save_dir}
}

function git_change_branch_ats()
{
  atsbranch=$1
  save_dir=`pwd`
  cd ${ats_source_dir}
  status_message "In ${ats_source_dir} checking out ATS branch ${atsbranch}"
  ${git_binary} checkout ${atsbranch}
  if [ $? -ne 0 ]; then
    error_message "Failed to update ATS at ${ats_source_dir} to branch ${atsbranch}"
    exit_now 30
  fi
  cd ${save_dir}
}

function git_submodule_clone()
{
  submodule_name=$1
  save_dir=`pwd`
  cd ${amanzi_source_dir}
  status_message "In ${amanzi_source_dir} checking out ${submodule_name}"
  ${git_binary} submodule update --init --remote ${submodule_name}
  if [ $? -ne 0 ]; then
    error_message "Failed to check out submodule ${submodule_name}"
    exit_now 30
  fi
  cd ${save_dir}
}

# ---------------------------------------------------------------------------- #
# MPI functions
# ---------------------------------------------------------------------------- #
function check_mpi_root
{
  if [ -z "${mpi_root_dir}" ]; then

    mpi_root_env="${MPIROOT} ${MPI_ROOT} ${MPIHOME} ${MPI_HOME} ${MPI_PREFIX}"
    for env_try in ${mpi_root_env}; do
      if [ -e "${env_try}" ]; then
        mpi_root_dir="${env_try}"
        status_message "Located MPI installation in ${mpi_root_dir}"
        break
      fi
    done

  else

    if [ ! -e "${mpi_root_dir}" ] ; then
      error_message "MPI root directory ${mpi_root_dir} does not exist"
      exit_now 30
    fi

  fi
}


# ---------------------------------------------------------------------------- #
# xSDK functions: check
# ---------------------------------------------------------------------------- #

function check_xsdk_root
{
  if [ -z "${xsdk_root_dir}" ]; then

    xsdk_root_env="${XSDKROOT} ${XSDK_ROOT} ${XSDKHOME} ${XSDK_HOME} ${XSDK_PREFIX}"
    for env_try in ${xsdk_root_env}; do
      if [ -e "${env_try}" ]; then
        xsdk_root_dir="${env_try}"
        status_message "Located xSDK installation in ${xsdk_root_dir}"
        break
      fi
    done

  else

    if [ ! -e "${xsdk_root_dir}" ] ; then
      error_message "xSDK root directory ${xsdk_root_dir} does not exist"
      status_message "Installing xSDK as a TPL"
    fi

  fi
}


# ---------------------------------------------------------------------------- #
# Spack functions: check
# ---------------------------------------------------------------------------- #

function check_Spack
{
  if [ ! -e "${Spack_binary}" ]; then
    error_message "Spack binary does not exist - Will try to locate..."

    if [ ! -e ${tpl_install_prefix}/spack/bin/spack ]; then
      error_message "Could not locate Spack - Downloading and installing as a TPL"
#     build_Spack=$TRUE
#     status_message "build_Spack: ${build_Spack}"

      pwd_save=`pwd`
      if [ ! -e ${tpl_install_prefix} ]; then
        mkdir_now ${tpl_install_prefix}
      fi
      cd ${tpl_install_prefix}
      git clone https://github.com/LLNL/spack.git
      

      if [ ${xsdk} == ${TRUE} ]; then
        cd ${tpl_install_prefix}/spack/bin
        #git checkout 9e95e83
        #git pull
        status_message "Installing xSDK..."
        ./spack install xsdk@0.3.0
      fi
      cd ${pwd_save}
    fi

      Spack_binary=${tpl_install_prefix}/spack/bin/spack
      status_message "Spack binary: ${Spack_binary}"

  else
    status_message "Spack binary: ${Spack_binary}"

  fi
}


# ---------------------------------------------------------------------------- #
# Compiler functions
# ---------------------------------------------------------------------------- #

function define_c_compiler
{
  if [ -z "${build_c_compiler}" ]; then
 
    status_message "Searching for a C compiler"

    # build a list to search
    c_search_list=
    if [ -n "${mpi_root_dir}" ]; then
      c_search_list="${mpi_root_dir}/bin/mpicc"
    fi
    
    if [ -n "${CC}" ]; then 
      c_search_list="${c_search_list} ${CC}"
    fi

    c_search_list="${c_search_list} ${known_c_compilers}"
    for c_try in ${c_search_list}; do
      status_message "Searching for ${c_try}"
      full_c_try=`which ${c_try} 2>/dev/null`
      if [ -e "${full_c_try}" ]; then
        build_c_compiler="${full_c_try}"
        break
      fi
    done

    if [ -z "${build_c_compiler}" ]; then
      error_message "Failed to locate a C compiler. Please use the --with-c-compiler option."
      exit_now 30
    fi

  fi

  status_message "Build with C compiler: ${build_c_compiler}"
  export CC=${build_c_compiler}
}

function define_cxx_compiler
{
  if [ -z "${build_cxx_compiler}" ]; then
 
    status_message "Searching for a C++ compiler"

    # build a list to search
    cxx_search_list=
    if [ -n "${mpi_root_dir}" ]; then
      cxx_search_list="${mpi_root_dir}/bin/mpiCC"
      cxx_search_list="${mpi_root_dir}/bin/mpicxx ${cxx_search_list}"
    fi
    
    if [ -n "${CXX}" ]; then 
      cxx_search_list="${cxx_search_list} ${CXX}"
    fi

    cxx_search_list="${cxx_search_list} ${known_cxx_compilers}"
    for cxx_try in ${cxx_search_list}; do
      status_message "Searching for ${cxx_try}"
      full_cxx_try=`which ${cxx_try} 2>/dev/null`
      if [ -e "${full_cxx_try}" ]; then
        build_cxx_compiler="${full_cxx_try}"
        break
      fi
    done

    if [ -z "${build_cxx_compiler}" ]; then
      error_message "Failed to locate a C++ compiler. Please use the --with-cxx-compiler option."
      exit_now 30
    fi

  fi

  status_message "Build with C++ compiler: ${build_cxx_compiler}"
  export CXX=${build_cxx_compiler}
}

function define_fort_compiler
{
   if [ -z "${build_fort_compiler}" ]; then
 
    status_message "Searching for a Fortran compiler"

    # build a list to search
    fort_search_list=
    if [ -n "${mpi_root_dir}" ]; then
      fort_search_list="${mpi_root_dir}/bin/mpif90"
    fi
    
    if [ -n "${FC}" ]; then 
      fort_search_list="${fort_search_list} ${FC}"
    fi

    fort_search_list="${fort_search_list} ${known_fortran_compilers}"
    for fort_try in ${fort_search_list}; do
      status_message "Searching for ${fort_try}"
      full_fort_try=`which ${fort_try} 2>/dev/null`
      if [ -e "${full_fort_try}" ]; then
        build_fort_compiler="${full_fort_try}"
        break
      fi
    done

    if [ -z "${build_fort_compiler}" ]; then
      error_message "Failed to locate a Fortran compiler. Please use the --with-fort-compiler option."
      exit_now 30
    fi

  fi

  status_message "Build with Fortran compiler: ${build_fort_compiler}"
  export FC=${build_fort_compiler}

}

function check_compilers
{
  define_c_compiler
  define_cxx_compiler
  define_fort_compiler

  status_message "Compiler Check complete"
}

version_compare_element()
{
  (( $1 == $2 )) && return 0 #equal
  (( $1 > $2 )) && return 1  #greater than
  (( $1 < $2 )) && return 2  #less than
  exit 1
}

version_compare()
{
  A=(${1//./ })
  B=(${2//./ })
  i=0
  while (( i < ${#A[@]} )) && (( i < ${#B[@]})); 
  do
    version_compare_element "${A[i]}" "${B[i]}"
    result=$?
    [[ $result =~ [12] ]] && return $result
    let i++
  done
  version_compare_element "${#A[i]}" "${#B[i]}"
  return $?
}


# ---------------------------------------------------------------------------- #
# Tools functions: checks
# ---------------------------------------------------------------------------- #

function check_tools
{
  # Check Git
  if [ ! -e "${git_binary}" ]; then
    error_message "Git binary does not exist"
    exit_now 10
  fi
  status_message "Git binary: ${git_binary}"

  # Check CURL 
  if [ ! -e "${curl_binary}" ]; then
    error_message "CURL binary does not exist"
    exit_now 10
  fi
  status_message "CURL binary: ${curl_binary}"

  # Check CMake and CTest
  if [ -z "${cmake_binary}" ]; then 
    status_message "CMake not defined. Will build"
    build_cmake ${tools_build_dir} ${tools_install_prefix} ${tools_download_dir} 
  fi
  if [ ! -e "${cmake_binary}" ]; then
    error_message "CMake binary does not exist. Will build."
    build_cmake ${tools_build_dir} ${tools_install_prefix} ${tools_download_dir} 
  else
    # CMake binary does exist, not check the version number
    ver_string=`${cmake_binary} --version | sed 's/cmake version[ ]*//'`
    status_message "Found CMake version:  ${ver_string}"
    ( version_compare "$ver_string" "$cmake_version" )
    result=$?
    if [ ${result} -eq 2 ]; then
      status_message "CMake version is less than required version. Will build CMake version 3.11.4"
      build_cmake ${tools_build_dir} ${tools_install_prefix} ${tools_download_dir} 
    fi
  fi
  if [ ! -e "${ctest_binary}" ]; then
    error_message "CTest binary does not exist. Will deactivate the test suite."
    test_suite=${FALSE}
    reg_tests=${FALSE}
  fi

  status_message "CMake binary: ${cmake_binary}"
  status_message "CTest binary: ${ctest_binary}"

  status_message "Tools check complete"
}


# ---------------------------------------------------------------------------- #
# Directory functions
# ---------------------------------------------------------------------------- #

function define_build_directories
{
  # The amanzi build directory
  if [ ! -e "${amanzi_build_dir}" ]; then
    mkdir_now "${amanzi_build_dir}"
  fi

  # The TPL build directory
  if [ ! -e "${tpl_build_dir}" ]; then
    mkdir_now "${tpl_build_dir}"
  fi

  # The TPL download directory
  if [ ! -e "${tpl_download_dir}" ]; then
    mkdir_now "${tpl_download_dir}"
  fi

  # Tools build directory
  if [ ! -e "${tools_build_dir}" ]; then
    mkdir_now "${tools_build_dir}"
  fi

  # Tools download directory
  if [ ! -e "${tools_download_dir}" ]; then
    mkdir_now "${tools_download_dir}"
  fi

  status_message "Build directories ready"
}

function define_install_directories
{
  # The prefix option overrides the other install choices
  if [ ! -z  "${prefix}" ] ; then
    status_message "Global prefix (${prefix}) set override Amanzi and TPL installations"
    amanzi_install_prefix=${prefix}
    tpl_install_prefix=${prefix}/tpls
  fi

  status_message "Amanzi installation: ${amanzi_install_prefix}"
  status_message "TPL installation: ${tpl_install_prefix}"
}    

function define_unstructured_dependencies
{
  if [ "${unstructured}" -eq "${FALSE}" ]; then
    eval "mstk_mesh=$FALSE"
    eval "moab_mesh=$FALSE"
  fi
}

function define_structured_dependencies
{
  if [ "${structured}" -eq "${TRUE}" ]; then
    eval "petsc=$TRUE"
    status_message "Enable package PETSc"
  fi
}


# ---------------------------------------------------------------------------- #
# Arch-specific functions
# ---------------------------------------------------------------------------- #
function define_nersc_options
{
  shared=$FALSE
  prefer_static=$TRUE
  exec_static=$TRUE
  prg_env="gnu"
    
  libsci_file=${tpl_build_src_dir}/include/trilinos-blas-libsci-${prg_env}.cmake
  arch_tpl_opts="-DAMANZI_ARCH:STRING=${amanzi_arch} \
                 -DMPI_EXEC:STRING=srun \
                 -DMPI_EXEC_NUMPROCS_FLAG:STRING=-n \
                 -DPREFER_STATIC_LIBRARIES:BOOL=${prefer_static} \
                 -DBUILD_STATIC_EXECUTABLES:BOOL=${exec_static} \
                 -DTrilinos_Build_Config_File:FILEPATH=${libsci_file}"
  
  arch_amanzi_opts="-DTESTS_REQUIRE_MPIEXEC:BOOL=${TRUE} \
                    -DTESTS_REQUIRE_FULLPATH:BOOL=${TRUE}"

  echo "Setting ARCH for: NERSC"
  echo "ARCH TPL OPTS = " ${arch_tpl_opts}
  echo "ARCH AMANZI OPTS = " ${arch_amanzi_opts}
}

function define_summit_options
{
  if [ "${kokkos_cuda}" -eq "${TRUE}" ]; then
    mpi_exec_args="-a 1 -c 1 -g 1" # 1 gpu per mpi rank, max 6 ranks
  elif [ "${kokkos_openmp}" -eq "${TRUE}" ]; then
    mpi_exec_args="-a 1 -c 7" # 7 cpus per mpi rank, max 6 ranks
  fi
  
  arch_tpl_opts="-DAMANZI_ARCH:STRING=${amanzi_arch} \
                 -DMPI_EXEC:STRING=jsrun \
                 -DMPI_EXEC_NUMPROCS_FLAG:STRING=-n \
                 -DMPI_EXEC_MAX_NUMPROCS:INT=6"
  arch_amanzi_opts="-DAMANZI_ARCH:STRING=${amanzi_arch} \
                 -DMPI_EXEC:STRING=jsrun \
                 -DMPI_EXEC_NUMPROCS_FLAG:STRING=-n \
                 -DMPI_EXEC_MAX_NUMPROCS:INT=6"
  echo "Setting ARCH for: Summit"
  echo "ARCH TPL OPTS = " ${arch_tpl_opts}
  echo "ARCH AMANZI OPTS = " ${arch_amanzi_opts}
}    


# ---------------------------------------------------------------------------- #
#  Main 
# ---------------------------------------------------------------------------- #

# Parse the command line arguments
array=( "$@" )
parse_argv "${array[@]}"
print_variable_values

# Define the TPL build source directory
tpl_build_src_dir=${amanzi_source_dir}/config/SuperBuild

# Set packages that depend on unstructured
define_unstructured_dependencies

# Set packages that depend on structured
define_structured_dependencies

# Define the install directories 
define_install_directories

# Create the build directories
define_build_directories

# Define the compilers
check_compilers

# Check the cmake, git and curl tools
check_tools

# Set extra options for building on nersc
if [ "${amanzi_arch}" = "Summit" ]; then
    define_summit_options
elif [ "${amanzi_arch}" = "NERSC" ]; then
    define_nersc_options
elif [ "${amanzi_arch}" != "" ]; then
    error_message "ARCH ${amanzi_arch} not supported -- valid are Summit and NERSC"
    exit_now 10
fi
    



# Print and exit if --print is set
if [ "${print_exit}" -eq "${TRUE}" ]; then
  print_variable_values
  exit_now 0
fi

# Change the branch
if [ ! -z "${amanzi_branch}" ]; then
  git_change_branch ${amanzi_branch}
fi

# ---------------------------------------------------------------------------- #
# Configure tools
# ---------------------------------------------------------------------------- #
if [ ! -n "${mpi_root_dir}" ]; then
  tools_build_src_dir=${amanzi_source_dir}/config/ToolsBuild

  cd ${tools_build_dir}
  ${cmake_binary} \
      -DCMAKE_C_FLAGS:STRING="${build_c_flags}" \
      -DCMAKE_CXX_FLAGS:STRING="${build_cxx_flags}" \
      -DCMAKE_Fortran_FLAGS:STRING="${build_fort_flags}" \
      -DCMAKE_EXE_LINKER_FLAGS:STRING="${build_link_flags}" \
      -DCMAKE_BUILD_TYPE:STRING=${tpls_build_type} \
      -DCMAKE_C_COMPILER:FILEPATH=${build_c_compiler} \
      -DCMAKE_CXX_COMPILER:FILEPATH=${build_cxx_compiler} \
      -DCMAKE_Fortran_COMPILER:FILEPATH=${build_fort_compiler} \
      -DTOOLS_MPI:STRING="${tools_mpi}" \
      -DTOOLS_INSTALL_PREFIX:STRING=${tools_install_prefix} \
      -DTOOLS_DOWNLOAD_DIR:FILEPATH=${tools_download_dir} \
      -DTOOLS_PARALLEL_JOBS:INT=${parallel_jobs} \
      ${tools_build_src_dir}

  if [ $? -ne 0 ]; then
    error_message "Failed to build tools"
    exit_now 30
  fi
  status_message "Tools configure complete"
 
  # Tools parameters
  # OpenMPI 3.x requires more slots to run Amanzi tests
  if [ "${tools_mpi}" = "openmpi" ]; then
    mpi_exec_args+="--oversubscribe"
  fi 
  
  # Tools install
  cd ${tools_build_dir}
  make install
  if [ $? -ne 0 ]; then
    error_message "Failed to install configure script"
    exit_now 30
  fi
      
  tools_config_file=${tools_install_prefix}/share/cmake/amanzi-tools-config.cmake
  mpi_root_dir=${tools_install_prefix}
  build_c_compiler=${tools_install_prefix}/bin/mpicc
  build_cxx_compiler=${tools_install_prefix}/bin/mpicxx
  build_fort_compiler=${tools_install_prefix}/bin/mpif90

 cd ${pwd_save}
      
  status_message "Tools build complete: new MPI root=${mpi_root_dir}"
  status_message "  new MPI root=${mpi_root_dir}"
  status_message "  new C compiler=${build_c_compiler}"
  status_message "For future Amanzi builds use ${tools_config_file}"
fi
  

# Check the MPI root value 
check_mpi_root

# Define the compilers
check_compilers


# ---------------------------------------------------------------------------- #
# Now build the TPLs if the config file is not defined
# ---------------------------------------------------------------------------- #
if [ -z "${tpl_config_file}" ]; then

  # Return to this directory once the TPL build is complete
  pwd_save=`pwd`

  # Check for Spack
  if [ "${Spack}" -eq "${TRUE}" ]; then
    status_message "Building with Spack"
    check_Spack
  fi

  # Are we using xSDK?  If so, make sure Spack in enabled
  if [ "${xsdk}" -eq "${TRUE}" ]; then 
    status_message "Building with xSDK"

    if [ "${Spack}" -eq "${FALSE}" ]; then
      status_message "Enabling Spack"
      Spack=$TRUE
      check_Spack
    fi
  fi

  if [ "${build_amanzi}" -eq "${FALSE}" ]; then
    status_message "Only building TPLs, stopping before building Amanzi itself"
  fi
 
 
  # Configure the TPL build
  cmd_configure="${cmake_binary} \
      -DCMAKE_C_FLAGS:STRING="${build_c_flags}" \
      -DCMAKE_CXX_FLAGS:STRING="${build_cxx_flags}" \
      -DCMAKE_Fortran_FLAGS:STRING="${build_fort_flags}" \
      -DCMAKE_EXE_LINKER_FLAGS:STRING="${build_link_flags}" \
      -DCMAKE_BUILD_TYPE:STRING=${tpls_build_type} \
      -DTRILINOS_BUILD_TYPE:STRING=${trilinos_build_type} \
      -DCMAKE_C_COMPILER:FILEPATH=${build_c_compiler} \
      -DCMAKE_CXX_COMPILER:FILEPATH=${build_cxx_compiler} \
      -DCMAKE_Fortran_COMPILER:FILEPATH=${build_fort_compiler} \
      -DMPI_PREFIX:STRING="${mpi_root_dir}" \
      -DTPL_INSTALL_PREFIX:STRING=${tpl_install_prefix} \
      -DENABLE_Structured:BOOL=${structured} \
      -DENABLE_Unstructured:BOOL=${unstructured} \
      -DENABLE_CCSE_TOOLS:BOOL=${ccse_tools} \
      -DCCSE_BL_SPACEDIM:INT=${spacedim} \
      -DENABLE_MOAB_Mesh:BOOL=${moab_mesh} \
      -DENABLE_MSTK_Mesh:BOOL=${mstk_mesh} \
      -DENABLE_NetCDF4:BOOL=${netcdf4} \
      -DENABLE_HYPRE:BOOL=${hypre} \
      -DENABLE_SUPERLU:BOOL=${superlu} \
      -DENABLE_PETSC:BOOL=${petsc} \
      -DENABLE_ALQUIMIA:BOOL=${alquimia} \
      -DENABLE_PFLOTRAN:BOOL=${pflotran} \
      -DENABLE_CRUNCHTOPE:BOOL=${crunchtope} \
      -DENABLE_Silo:BOOL=${silo} \
      -DENABLE_SPACK:BOOL=${Spack} \
      -DENABLE_EPETRA:BOOL=${epetra} \
      -DENABLE_KOKKOS:BOOL=${kokkos} \
      -DENABLE_KOKKOS_CUDA:BOOL=${kokkos_cuda} \
      -DENABLE_KOKKOS_OPENMP:BOOL=${kokkos_openmp} \
      -DSPACK_BINARY:STRING=${Spack_binary} \
      -DBUILD_SPACK:BOOL=${build_Spack} \
      -DENABLE_XSDK:BOOL=${xsdk} \
      -DBUILD_SHARED_LIBS:BOOL=${shared} \
      -DTPL_DOWNLOAD_DIR:FILEPATH=${tpl_download_dir} \
      ${arch_tpl_opts} \
      ${tpl_build_src_dir}"

  # TPL build will create this configuration file
  tpl_config_file=${tpl_install_prefix}/share/cmake/amanzi-tpl-config.cmake
  
  # Echo or execute configure command
  if [ ${dry_run} == "${TRUE}" ] ; then
    status_message "TPL configure command: \n"
    echo ${cmd_configure}
    echo ""
  else
    cd ${tpl_build_dir}
    ${cmd_configure}

    if [ $? -ne 0 ]; then
      error_message "Failed to configure TPL build"
      exit_now 30
    fi
  fi

  if [ ${dry_run} == "${FALSE}" ]; then
    status_message "TPL configure complete"
  
    # TPL make 
    cd ${tpl_build_dir}
    make -j ${parallel_jobs}
    if [ $? -ne 0 ]; then
      error_message "Failed to build TPLs"
      exit_now 30
    fi
  
    # TPL Install
    cd ${tpl_build_dir}
    make install
    if [ $? -ne 0 ]; then
      error_message "Failed to install configure script"
      exit_now 30
    fi
      
    cd ${pwd_save}
      
    status_message "TPL build complete"
    status_message "For future Amanzi builds use ${tpl_config_file}"

  else
    status_message "To execute this TPL build remove the --dry_run option."
  fi
  
else 
  status_message "Checking configuration file ${tpl_config_file}"

  if [ ! -e "${tpl_config_file}" ]; then
    error_message "Configure file ${tpl_config_file} does not exist"
    exit_now 30
  fi

  if [ ! -r "${tpl_config_file}" ]; then
    error_message "Configure file ${tpl_config_file} is not readable"
    exit_now 30
  fi

  if [ ! -f "${tpl_config_file}" ]; then
    error_message "Configure file ${tpl_config_file} is not a regular file"
    exit_now 30
  fi

fi

if [ "${build_amanzi}" -eq "${FALSE}" ]; then
    exit_now 0
fi

status_message "Build Amanzi with configure file ${tpl_config_file}"

# Clone any submodules (ATS)
if [ "${ats_physics}" -eq "${TRUE}" ]; then
    if [ "${amanzi_physics}" -eq "${TRUE}" ]; then
        error_message "amanzi_physics and ats_physics are currently incompatible -- enable only one."
        exit_now 30
    fi

    git_submodule_clone "src/physics/ats"
    if [ ! -z "${ats_branch}" ]; then
        git_change_branch_ats ${ats_branch}
    fi
fi


# Configure the Amanzi build
cmd_configure="${cmake_binary} \
    -C${tpl_config_file} \
    -DCMAKE_C_FLAGS:STRING="${build_c_flags}" \
    -DCMAKE_CXX_FLAGS:STRING="${build_cxx_flags}" \
    -DCMAKE_Fortran_FLAGS:STRING="${build_fort_flags}" \
    -DCMAKE_EXE_LINKER_FLAGS:STRING="${build_link_flags}" \
    -DCMAKE_INSTALL_PREFIX:STRING=${amanzi_install_prefix} \
    -DCMAKE_BUILD_TYPE:STRING=${build_type} \
    -DENABLE_Structured:BOOL=${structured} \
    -DENABLE_Unstructured:BOOL=${unstructured} \
    -DENABLE_MOAB_Mesh:BOOL=${moab_mesh} \
    -DENABLE_MSTK_Mesh:BOOL=${mstk_mesh} \
    -DENABLE_SUPERLU:BOOL=${superlu} \
    -DENABLE_HYPRE:BOOL=${hypre} \
    -DENABLE_PETSC:BOOL=${petsc} \
    -DENABLE_ALQUIMIA:BOOL=${alquimia} \
    -DENABLE_PFLOTRAN:BOOL=${pflotran} \
    -DENABLE_CRUNCHTOPE:BOOL=${crunchtope} \
    -DENABLE_Silo:BOOL=${silo} \
    -DENABLE_EPETRA:BOOL=${epetra} \
    -DENABLE_KOKKOS:BOOL=${kokkos} \
    -DENABLE_KOKKOS_CUDA:BOOL=${kokkos_cuda} \
    -DENABLE_KOKKOS_OPENMP:BOOL=${kokkos_openmp} \
    -DENABLE_AmanziPhysicsModule:BOOL=${amanzi_physics} \
    -DENABLE_ATSPhysicsModule:BOOL=${ats_physics} \
    -DBUILD_SHARED_LIBS:BOOL=${shared} \
    -DCCSE_BL_SPACEDIM:INT=${spacedim} \
    -DENABLE_Regression_Tests:BOOL=${reg_tests} \
    -DMPI_EXEC_GLOBAL_ARGS:STRING=${mpi_exec_args} \
    ${arch_amanzi_opts} \
    ${amanzi_source_dir}"

# Echo or execute configure command
if [ ${dry_run} == "${TRUE}" ] ; then
   status_message "Amanzi configure command: \n"
   echo ${cmd_configure}
   echo ""
else
   cd ${amanzi_build_dir}
   ${cmd_configure}
fi

if [ ${dry_run} == "${FALSE}" ]; then
  if [ $? -ne 0 ]; then
    error_message "Failed to configure Amanzi"
    exit_now 50
  fi
  status_message "Amanzi configure complete"

  # Amanzi Build
  if [ "${build_stage_1}" == "${TRUE}" ]; then
    # we could make this fancier, but operators is roughly mid-way in the build in terms of dependencies.
    status_message "Amanzi build-stage 1: "
    cd ${amanzi_build_dir}/src/whetstone
    status_message "  - building in $(pwd)"
    make -j ${parallel_jobs}
    cd ${amanzi_build_dir}/src/operators
    status_message "  - building in $(pwd)"
    make -j ${parallel_jobs}
    cd ${amanzi_build_dir}
    if [ $? -ne 0 ]; then
        error_message "Amanzi build-stage 1: FAILED"
        exit_now 50
    else
        status_message "Amanzi build-stage 1: complete"
        status_message "Bootstrap Incomplete: run with --enable-build_stage_2 to complete the build."
        exit_now 0
    fi
  elif [ "${build_stage_2}" == "${TRUE}" ]; then
    status_message "Amanzi build-stage 2: "
    status_message "  - building in $(pwd)"
    make -j ${parallel_jobs}
    if [ $? -ne 0 ]; then
        error_message "Amanzi build-stage 2: FAILED"
        exit_now 50
    else
        status_message "Amanzi build-stage 2: complete"
    fi
  else
    make -j ${parallel_jobs}
    if [ $? -ne 0 ]; then
        error_message "Failed to build Amanzi"
        exit_now 50
    fi
    status_message "Amanzi build complete"
  fi

  # Amanzi Test Suite
  if [ "${test_suite}" -eq "${TRUE}" ]; then
    status_message "Run Amanzi test suite"
    ${ctest_binary} --output-on-failure --output-log amanzi-test-results.log #-j ${parallel_jobs}
    status_message "Test results in amanzi-test-results.log"
    if [ $? -ne 0 ]; then
      error_message "Amanzi test suite failed"
      exit_now 30
    fi
  fi

  # Amanzi Install
  if [ "${build_stage_1}" -ne "${TRUE}" ]; then
    make install
  fi
  if [ $? -ne 0 ]; then
    error_message "Failed to install Amanzi"
    exit_now 50
  fi
  status_message "Amanzi install complete"

else
  status_message "To execute this Amanzi build remove the --dry_run option."
fi

status_message "Build Amanzi with configure file ${tpl_config_file}"


# Build User Guide
if [ ${dry_run} == "${FALSE}" ]; then
  if [ "${build_user_guide}" -eq "${FALSE}" ]; then
    exit_now 0
  fi

  cd ${amanzi_source_dir}/doc/user_guide

  export AMANZI_INSTALL_DIR=${amanzi_install_prefix}
  export AMANZI_MPI_EXEC=${mpi_root_dir}/bin/mpiexec
  export AMANZI_MPI_NP=${parallel_jobs}
  export PYTHONPATH=${tpl_install_prefix}/lib:$(eval echo \$\{PYTHONPATH\})

  make html AMANZI_SECTION_OPTIONS="--full-guide"

  status_message "Build User Guide, see ${amanzi_source_dir}/doc/user_guide/_build/html/index.html"
fi

status_message "Bootstrap complete"

exit_now 0

