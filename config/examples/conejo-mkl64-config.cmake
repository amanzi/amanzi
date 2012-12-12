# ############################################################################ #
#
# Third Party Library (TPL Build)
#  CMake configuration file
#
# Machine: Conejo (LANL)
# OS: Linux Red Hat EL5 x86_64
# Compiler: Intel 11.1 with MKL 10.3 libraries
# MPI: OpenMPI 1.4.3
#
# Usage:
#   (1) Load the default Intel, OpenMPI and  MKL modules
#        module load intel openmpi-intel mkl/10.3
#   (2) Configure 
#       cmake -C <root path>/conejo-config.cmake <root directory amanzi>/amanzi/config/SuperBuild
# 
# ############################################################################ #

# --- Machine specific directories
set(ASCEM_PROJECT_DIR "/usr/projects/ascem" CACHE PATH "ASCEM project directory")
set(MPI_ROOT "$ENV{MPI_ROOT}" CACHE PATH "MPI installation location")
set(MKL_ROOT "$ENV{MKLROOT}" CACHE PATH "MKL installation location")
set(MKL_ARCH "intel64" CACHE STRING "MKL arch type")

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
# Vendor string Intel10_64lp 64-bit MKL libraries. ONLY vendor string that found correct libraries.
set(BLA_VENDOR "Intel10_64lp" CACHE STRING "Search for 64-bit MKL libraries")

# --- Mesh Capabilities
# Structured mesh code does not compile with Intel Fortran 
set(ENABLE_Structured FALSE CACHE BOOL "Flag for structured mesh capability")

# --- Solver Capabilities
#
set(ENABLE_HYPRE TRUE CACHE BOOL "Enable the HYPRE preconditioner package")

# --- TPL Installation location
set(TPL_INSTALL_PREFIX "${ASCEM_PROJECT_DIR}/tpls/installs/openmpi-intel-mkl64" CACHE PATH "ASCEM TPL installation location")
