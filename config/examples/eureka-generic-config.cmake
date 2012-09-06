# ############################################################################ #
#
# Third Party Library (TPL Build)
#  CMake configuration file
#
# Machine: Eureka (NSCEE)
# OS: Linux Red Hat EL5 x86_64
# Compiler: GCC 4.6.3
# MPI: PMPI (HP-MPI)
#
# Usage:
#   (1) Load the modules: cmake/cmake-2.8.8 PMPI/mdulefile gnu/gcc-4.6.3 
#        module load intel openmpi-intel atlas
#   (2) Configure 
#       cmake -C <root path>/eureka-atlas-config.cmake <root directory amanzi>/amanzi/config/SuperBuild
# 
# ############################################################################ #

# --- Machine specific directories
set(ASCEM_PROJECT_DIR "/home/ASCEM" CACHE PATH "ASCEM project directory")
set(MPI_ROOT "$ENV{MPI_ROOT}" CACHE PATH "MPI installation location")

# --- Set the build type Release == minimal -O3 optimization
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "CMake build type") 

# --- Set the compilers
set(CMAKE_C_COMPILER "${MPI_ROOT}/bin/mpicc" CACHE FILEPATH "MPI C compiler wrapper" FORCE)
set(CMAKE_CXX_COMPILER "${MPI_ROOT}/bin/mpiCC" CACHE FILEPATH "MPI C++ compiler wrapper" FORCE)
set(CMAKE_Fortran_COMPILER "${MPI_ROOT}/bin/mpif90" CACHE FILEPATH "MPI Fortran compiler wrapper" FORCE)

# --- Set test suite behavior
# Like many clusters, parallel runs are not allowed on login/compiler nodes.
# So we explicitly define the MPI executable and prevent the test suite from
# searching and finding the incorrect binary.
set(MPI_EXEC "${MPI_ROOT}/bin/mpirun" CACHE FILEPATH "MPI execute command")
set(MPI_EXEC_NUMPROCS_FLAG "-np" CACHE FILEPATH "MPI number of processes flag")

# --- LAPACK/BLAS Definitions
# Use the system lapack blas. Need static libs because the shared ones are broken.
# See symbol look up errors at run time. 
set(BLA_VENDOR "Generic" CACHE STRING "BLAS/LAPACK vendor search type")
set(BLA_STATIC TRUE CACHE BOOL "Search for static libraries")

# --- Solver Capabilities
set(ENABLE_HYPRE TRUE CACHE BOOL "Activate HYPRE APIs in Trilinos build")

# --- Mesh Capabilities

# --- TPL Installation location
set(TPL_INSTALL_PREFIX "${ASCEM_PROJECT_DIR}//tpls/installs/PMPI-gnu46/generic" CACHE PATH "ASCEM TPL installation location")
