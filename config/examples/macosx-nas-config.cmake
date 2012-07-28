# ############################################################################ #
#
# Third Party Library (TPL Build)
#  CMake configuration file
#
# Machine: Laptop
# OS: MacOSX 10.6.8
# Compiler: GCC 4.6 installed through MacPorts
# MPI: OpenMPI 1.4.4 installed from source
#
# Usage:
#   Configure 
#       cmake -C <root path>/macosx-nas-config.cmake <root directory amanzi>/amanzi/config/SuperBuild
# 
# ############################################################################ #

# --- Machine specific directories
set(MPI_ROOT "$ENV{MPI_ROOT}" CACHE PATH "MPI installation location")

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

# --- LAPACK/BLAS Definitions
set(ENABLE_BLA_Search TRUE CACHE BOOL "Activate the CMake BLAS/LAPACK search")
set(BLA_VENDOR "NAS" CACHE STRING "Search for the Apple vecLib libraries")


# --- Solver Capabilities
set(ENABLE_HYPRE TRUE CACHE BOOL "Flag to activate HYPRE build")

# --- Mesh Capabilities
# Structured mesh code does not compile with Intel Fortran 
set(ENABLE_Structured FALSE CACHE BOOL "Flag for structured mesh capability")

# --- TPL Installation location
set(TPL_INSTALL_PREFIX "$ENV{HOME}/projects/ascem/tpls/openmpi-gcc46-nas" CACHE PATH "ASCEM TPL installation location")
