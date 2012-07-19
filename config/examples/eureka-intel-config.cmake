# ############################################################################ #
#
# Third Party Library (TPL Build)
#  CMake configuration file
#
# Machine: Eureka (NSCEE)
# OS: Linux Red Hat EL5 x86_64
# Compiler: Intel 12.X with MKL (arch=intel64, 4-byte ints)
# MPI: MPICH
#
# Usage:
#   (1) Load the modules: cmake/cmake-2.8.8 intel/inteli-12-impi 
#   (2) Configure 
#       cmake -C <root path>/eureka-intel-config.cmake <root directory amanzi>/amanzi/config/SuperBuild
#
# To build Amanzi against this installation you will need to add the MPI library path
# under the INTEL MPI to LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=${I_MPI_ROOT}/intel64/lib:${LD_LIBRARY_PATH}
# or
# setenv LD_LIBRARY_PATH ${I_MPI_ROOT}/intel64/lib:${LD_LIBRARY_PATH}
# 
# ############################################################################ #

# --- Machine specific directories
set(ASCEM_HOME_DIR "/home/ASCEM" CACHE PATH "ASCEM home (project) directory")
set(INTEL_ARCH intel64 CACHE STRING "Intel 64-bit arch")
set(MPI_ROOT "$ENV{I_MPI_ROOT}/${INTEL_ARCH}" CACHE PATH "MPI installation location")
set(MKL_ROOT "$ENV{MKLROOT}" CACHE PATH "MKL installation location")

# --- Set the build type Release == minimal -O3 optimization
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "CMake build type") 

# --- Set the compilers
set(CMAKE_C_COMPILER "${MPI_ROOT}/bin/mpiicc" CACHE FILEPATH "MPI C compiler wrapper" FORCE)
set(CMAKE_CXX_COMPILER "${MPI_ROOT}/bin/mpiicpc" CACHE FILEPATH "MPI C++ compiler wrapper" FORCE)
set(CMAKE_Fortran_COMPILER "${MPI_ROOT}/bin/mpiifort" CACHE FILEPATH "MPI Fortran compiler wrapper" FORCE)

# --- Set test suite behavior
# Like many clusters, parallel runs are not allowed on login/compiler nodes.
# So we explicitly define the MPI executable and prevent the test suite from
# searching and finding the incorrect binary.
set(MPI_EXEC "${MPI_ROOT}/bin/mpirun" CACHE FILEPATH "MPI execute command")
set(MPI_EXEC_NUMPROCS_FLAG "-n" CACHE FILEPATH "MPI number of processes flag")

# --- LAPACK/BLAS Definitions
# I used the online tool http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
# to define the link options. I selected the options: intel64(arch), 32bit ints, multi-threaded libraries
# and dynamic linking against a single file. The static linking broke with undefined links
# that *should be found in libmkl_intel_thread.
set(BLAS_LIBRARIES
   "-L${MKL_ROOT}/lib/${INTEL_ARCH};-lmkl_rt;-lpthread;-lm" 
    CACHE STRING "MKL dynamic link BLAS libraries")

set(LAPACK_LIBRARIES
   "-L${MKL_ROOT}/lib/${INTEL_ARCH};-lmkl_rt;-lpthread;-lm" 
    CACHE STRING "MKL dynamic link LAPACK libraries")

# --- Solver Capabilities
set(ENABLE_HYPRE TRUE CACHE BOOL "Activate HYPRE APIs in Trilinos build")

# --- Mesh Capabilities
# CCSE can not compile with GCC 4.1.2
set(Enable_Structured FALSE CACHE BOOL "Deactivate structured mesh capabilities")

# --- Compiler flag options
# The MPICH_IGNORE_CXX_SEEK is required to avoid an MPICH SEEK_SET name clash with stdio.h
# The other flags are from the -fast documentation. Use all flags from -fast that do not require 
# static linking or ip0, which requires a special AR definition and the SuperBuild is not currently
# designed to handle that.
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build release (optimized) libraries and binaries")
set(CMAKE_C_FLAGS "-DMPICH_IGNORE_CXX_SEEK -O3 -no-prec-div -xHost" CACHE STRING "Intel -fast compiler option")
set(CMAKE_CXX_FLAGS "-DMPICH_IGNORE_CXX_SEEK -O3 -no-prec-div -xHost" CACHE STRING "Intel -fast compiler option")
set(CMAKE_Fortran_FLAGS "-O3 -no-prec-div -xHost" CACHE STRING "Intel -fast compiler option")

#set(CMAKE_AR "/shared/local/opt/intel/composer_xe_2011_sp1.6.233/bin/intel64/xiar" CACHE FILEPATH "Intel AR")

# --- TPL Installation location
set(TPL_INSTALL_PREFIX "${ASCEM_HOME_DIR}/tpls/installs/${INTEL_ARCH}-optimized"  CACHE PATH "ASCEM TPL installation location")
