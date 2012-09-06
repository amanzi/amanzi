# ############################################################################ #
#
# Third Party Library (TPL Build)
#  CMake configuration file
#
# Machine: Eureka
# OS: Linux Red Hat EL6 x86_64
# Compiler: Intel 12 with MKL 10.3 libraries
# MPI: 
#
# Usage:
#   (1) Load modules
#        module load cmake intel/intel-12 PMPI/modulefile
#   (2) Configure 
#       cmake -C <root path>/eureka-mkl64-config.cmake <root directory amanzi>/amanzi/config/SuperBuild
# 
# ############################################################################ #

# --- Machine specific directories
set(ASCEM_PROJECT_DIR "/home/ASCEM" CACHE PATH "ASCEM project directory")
set(MPI_ROOT "$ENV{MPI_ROOT}" CACHE PATH "MPI installation location")
set(MKL_ROOT "$ENV{MKLROOT}" CACHE PATH "MKL installation location")

# --- Set the build type Release == minimal -O3 optimization
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "CMake build type") 

# --- Set the compilers
set(CMAKE_C_COMPILER "icc" CACHE FILEPATH "MPI C compiler wrapper" FORCE)
set(CMAKE_CXX_COMPILER "icpc" CACHE FILEPATH "MPI C++ compiler wrapper" FORCE)
set(CMAKE_Fortran_COMPILER "ifort" CACHE FILEPATH "MPI Fortran compiler wrapper" FORCE)

# --- Set test suite behavior
# Like many clusters, parallel runs are not allowed on login/compiler nodes.
# So we explicitly define the MPI executable and prevent the test suite from
# searching and finding the incorrect binary.
set(MPI_INSTALL_PREFIX ${MPI_ROOT} CACHE PATH "MPI install location")
set(MPI_EXEC "${MPI_ROOT}/bin/mpirun" CACHE FILEPATH "MPI execute command")

# --- Download behavior

# --- LAPACK/BLAS Definitions
# Vendor string Intel10_64lp 64-bit MKL libraries. Will work if MKL_ROOT/lib is in
# LD_LIBRARY_PATH
set(BLA_VENDOR "Intel10_64lp" CACHE STRING "Search for 64-bit MKL libraries")

# --- Mesh Capabilities
# Structured mesh code does not compile with Intel Fortran 
set(ENABLE_Structured FALSE CACHE BOOL "Flag for structured mesh capability")

# --- Solver Capabilities
#
set(ENABLE_HYPRE TRUE CACHE BOOL "Enable the HYPRE preconditioner package")

# --- TPL Installation location
set(TPL_INSTALL_PREFIX "$ENV{HOME}/tpls/installs/PMPI-intel12/mkl64" CACHE PATH "ASCEM TPL installation location")
