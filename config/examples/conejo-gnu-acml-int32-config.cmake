# ############################################################################ #
#
# Third Party Library (TPL Build)
#  CMake configuration file
#
# Machine: Conejo (LANL)
# OS: Linux Red Hat EL5 x86_64
# Compiler: GCC 4.5.2 with ACML 32 bit ints
# MPI: OpenMPI 1.4.3
#
# Usage:
#   (1) Load the default GNU, OpenMPI and ACML modules
#        module load gcc openmpi-gcc acml-gcc
#   (2) Configure 
#       cmake -C <root path>/conejo-gnu-acml-int32-config.cmake <root directory amanzi>/amanzi/config/SuperBuild
# 
# ############################################################################ #

# --- Machine specific directories
set(ASCEM_PROJECT_DIR "/usr/projects/ascem" CACHE PATH "ASCEM project directory")
set(MPI_ROOT "$ENV{MPI_ROOT}" CACHE PATH "MPI installation location")
#set(AMCL_ROOT "/opt/ACML/acml-4.3.0/ifort64/lib" CACHE PATH "ACML installation location")

# --- Set the build type Release == minimal -O3 optimization
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "CMake build type") 

# --- Set the compilers
set(CMAKE_C_COMPILER "${MPI_ROOT}/bin/mpicc" CACHE FILEPATH "MPI C compiler wrapper" FORCE)
set(CMAKE_CXX_COMPILER "${MPI_ROOT}/bin/mpicxx" CACHE FILEPATH "MPI C++ compiler wrapper" FORCE)
set(CMAKE_Fortran_COMPILER "${MPI_ROOT}/bin/mpif90" CACHE FILEPATH "MPI Fortran compiler wrapper" FORCE)

# --- Set test suite behavior
# Like many clusters, parallel runs are not allowed on login/compiler nodes.
# So we explicitly define the MPI executable and prevent the test suite from
# searching and finding the incorrect binary.
set(MPI_EXEC "${MPI_ROOT}/bin/mpirun" CACHE FILEPATH "MPI execute command")

# --- Download behavior
# Conejo does not allow external downloads
set(DISABLE_EXTERNAL_DOWNLOAD TRUE CACHE BOOL "Disable external web site downloads")
set(TPL_DOWNLOAD_DIR "${ASCEM_PROJECT_DIR}/tpls/source_files" CACHE PATH "Location of the TPL distribution files")

# --- LAPACK/BLAS Definitions
set(ENABLE_BLA_Search TRUE CACHE BOOL "Activate the CMake BLAS/LAPACK search")
set(BLA_VENDOR "ACML" CACHE STRING "Set CMake vendor search to ACML")


# --- Solver Capabilities
set(ENABLE_HYPRE TRUE CACHE BOOL "Flag to activate HYPRE build")

# --- Mesh Capabilities
# Structured mesh code does not compile with Intel Fortran 
set(ENABLE_Structured FALSE CACHE BOOL "Flag for structured mesh capability")

# --- TPL Installation location
#set(TPL_INSTALL_PREFIX "${ASCEM_PROJECT_DIR}/tpls/installs/intel/acml/state-dev" CACHE PATH "ASCEM TPL installation location")
