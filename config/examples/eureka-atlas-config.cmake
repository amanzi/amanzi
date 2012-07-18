# ############################################################################ #
#
# Third Party Library (TPL Build)
#  CMake configuration file
#
# Machine: Eureka (NSCEE)
# OS: Linux Red Hat EL5 x86_64
# Compiler: GCC 4.1.2
# MPI: OpenMPI 1.4.?
#
# Usage:
#   (1) Load the modules: cmake/cmake-2.8.8 PMPI/mdulefile mpi/openmpi-interconnects-gnu 
#        module load intel openmpi-intel atlas
#   (2) Configure 
#       cmake -C <root path>/eureka-atlas-config.cmake <root directory amanzi>/amanzi/config/SuperBuild
# 
# ############################################################################ #

# --- Machine specific directories
set(MPI_ROOT "/opt/platform_mpi" CACHE PATH "MPI installation location")
set(ATLAS_ROOT "/usr/lib64/atlas" CACHE PATH "ATLAS installation location")

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
set(MPI_EXEC_NUMPROCS_FLAG "-np" CACHE FILEPATH "MPI number of processes flag")

# --- LAPACK/BLAS Definitions
# The LANL support recommends the following link flags when linking against ATLAS
# -L/usr/lib64 -lgfortran -L${ATLAS_ROOT}/lib -llapack -lf77blas -lcblas -latlas
# I used the same flags here because this machine has a similar OS.
# Since the gfortran library is shared, I left the -L,-l flags in the library string.
# The remaining libraries are all static; I used full paths for those library names. 
# Need the gfotran library to link to the BLAS libraries.
# We use the LAPACK built with ATLAS, thus we need the cblas library
# to link correctly to lapack. The ';' is the list delimiter in CMake. 
set(BLAS_LIBRARIES
    "-L/usr/lib64;-lgfortran;${ATLAS_ROOT}/libf77blas.a;${ATLAS_ROOT}/libatlas.a"
    CACHE STRING "ATLAS BLAS libraries")
set(LAPACK_LIBRARIES
     "-L/usr/lib64;-lgfortran;${ATLAS_ROOT}/liblapack.a;${ATLAS_ROOT}/libcblas.a;${ATLAS_ROOT}/libf77blas.a;${ATLAS_ROOT}/libatlas.a" 
   CACHE STRING "ATLAS LAPACK libraries")

# --- Solver Capabilities
set(ENABLE_HYPRE TRUE CACHE BOOL "Activate HYPRE APIs in Trilinos build")

# --- Mesh Capabilities
# CCSE can not compile with GCC 4.1.2
set(Enable_Structured FALSE CACHE BOOL "Deactivate structured mesh capabilities")

# --- TPL Installation location
set(TPL_INSTALL_PREFIX "$ENV{HOME}/ascem/install/tpls" CACHE PATH "ASCEM TPL installation location")
