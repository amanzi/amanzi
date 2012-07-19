# ############################################################################ #
#
# Third Party Library (TPL Build)
#  CMake configuration file
#
# Machine: Conejo (LANL)
# OS: Linux Red Hat EL5 x86_64
# Compiler: Intel 11.1
# MPI: OpenMPI 1.4.3
#
# Usage:
#   (1) Load the default Intel, OpenMPI and ATLAS modules
#        module load intel openmpi-intel atlas
#   (2) Configure 
#       cmake -C <root path>/conejo-config.cmake <root directory amanzi>/amanzi/config/SuperBuild
# 
# ############################################################################ #

# --- Machine specific directories
set(ASCEM_PROJECT_DIR "/usr/projects/ascem" CACHE PATH "ASCEM project directory")
set(MPI_ROOT "$ENV{MPI_ROOT}" CACHE PATH "MPI installation location")
set(ATLAS_ROOT "$ENV{ATLAS_ROOT}" CACHE PATH "ATLAS installation location")

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
# The LANL support recommends the following link flags when linking against ATLAS
# -L/usr/lib64 -lgfortran -L${ATLAS_ROOT}/lib -llapack -lf77blas -lcblas -latlas
# Since the gfortran library is shared, I left the -L,-l flags in the library string.
# The remaining libraries are all static, so I used full paths there. 
# Need the gfotran library to link to the BLAS libraries.
# We use the LAPACK built with ATLAS, thus we need the cblas library
# to link correctly to lapack. The ';' is the list delimiter in CMake. 
set(BLAS_LIBRARIES
    "-L/usr/lib64;-lgfortran;${ATLAS_ROOT}/lib/libf77blas.a;${ATLAS_ROOT}/lib/libatlas.a"
    CACHE STRING "ATLAS BLAS libraries")
set(LAPACK_LIBRARIES
     "-L/usr/lib64;-lgfortran;${ATLAS_ROOT}/lib/liblapack.a;${ATLAS_ROOT}/lib/libcblas.a;${ATLAS_ROOT}/lib/libf77blas.a;${ATLAS_ROOT}/lib/libatlas.a" 
   CACHE STRING "ATLAS LAPACK libraries")

# --- Mesh Capabilities
# Structured mesh code does not compile with Intel Fortran 
set(ENABLE_Structured FALSE CACHE BOOL "Flag for structured mesh capability")

# --- TPL Installation location
set(TPL_INSTALL_PREFIX "${ASCEM_PROJECT_DIR}/tpls/installs/openmpi-intel" CACHE PATH "ASCEM TPL installation location")
