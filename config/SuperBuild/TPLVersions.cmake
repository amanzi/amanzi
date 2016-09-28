#  -*- mode: cmake -*-

#
# TPLVersions
#
#    Define the versions, approved download locations for each TPL
#

#
# TPL: Amanzi Collection of TPLs
#
#   Define a "version number" for the collection of TPLs listed here.
#   It's not clear this is the best way to include this information, 
#   but it's a reasonable place to start.
#   
#   Upgrade History:
#
#   0.90.6       - first version reference used in installations
#   0.90.7       - updated MSTK to version 2.01
#                - added NETCDF - fortran version 4.2 (interface library)
#   0.90.8       - added Xerces-C++ version 3.1.1 (XML support)
#   0.90.9       - update MSTK to version 2.02
#
#   0.91.0       - added PFLOTRAN and Alquimia (updates from Jeff Johnson's work on state-branch)
#   0.91.1       - updated CCSE to version 1.1.7 (and added patch)
#   0.91.2       - turned on patch for IfPack support of noncontiguous global ids with HYPRE
#   0.91.3       - updated Xerces-C++ build/find to use OSX CoreServices framework
#   0.91.4       - updated CCSE to version 1.1.8
#   0.91.5       - updated CCSE to version 1.2.1
#   0.91.6       - updated CCSE to version 1.2.3
#   0.91.7       - updated MSTK to version 2.10rc3
#   0.91.8       - updated Alquimia to licensed version 0.1
#   0.91.9       - updated CCSE to version 1.2.4
#   0.91.10      - updated MSTK to version 2.10rc5
#   0.91.11      - updated PFlotran to commit 1afe88d.
#   0.91.12      - updated MSTK to version 2.10
#   0.91.12a     - updated METIS to 5.1.0, ParMetis to 4.0.3a, SuperLU to 4.3, SuperLUDist to 3.3 and PETSc to 3.4.3
#   0.91.12b     - updated Trilinos to 10.6.1 and MSTK to 2.11rc2
#   0.91.12c     - updated CCSE to version 1.2.5
#   0.91.13      - updated MSTK to version 2.11rc3
#   0.91.14      - updated MSTK to version 2.11rc4 (fixes memory leaks)
#   0.91.15      - updated MSTK to version 2.11rc5 (fixes memory leaks)
#
#   0.92.0       - Merge lib updates through 0.91.15
#   0.92.1       - update MSTK to version 2.12 (fixes debug version linking)
#   0.92.2       - update CCSE to version 1.2.7 (adds f90 utility for plotting)
#   0.92.3       - update CURL to version 7.37.0 (builds correctly on Mac OS X 10.9)
#   0.92.4       - update NetCDF to version 4.3.2 (builds correctly on Mac OS X 10.9)
#   0.92.5       - Patched Alquimia to build properly with GFortran 4.9.x
#   0.92.6       - update Boost to version 1.56.0 
#   0.92.7       - update CCSE to version 1.2.8 
#   0.92.8       - update ExodussII 5.22 -> 6.06
#   0.92.9       - update MSTK to v 2.21 (incompatible -DWITH_MSTK_2_21rc1_OR_NEWER=TRUE)
#   0.92.10      - update MSTK to v 2.22rc1
#   0.92.11      - update MSTK to v 2.22rc3 (fixed parallel mesh partitioning bug)
#   0.92.12      - update PETSc to 3.5.2, Alquimia to 0.2, and PFlotran to commit 611092f80ddb.
#   0.92.13      - update MSTK to v2.22, includes installation of mesh utilities
#   0.92.14      - update Hypre to v2.10.0b (and added patch for to ensure tol>0)
#   0.92.15      - updated Alquimia to v0.2 (backward compatible)
#   0.92.16      - update CCSE to version 1.3.0 
#   0.92.17      - update MSTK to version 2.23 (adds element set capabilities)
#   0.92.18      - update Boost to version 1.58.0 
#   0.92.19      - update CCSE to version 1.3.2
#   0.92.20      - update CCSE to version 1.3.4 (fix issue with fsnapshot)
#   0.92.21      - update Alquimia to version 0.3.1 (CrunchFlow integration)
#   0.92.22      - added optional Silo package
#   0.92.23      - Patched ASCEM-IO to allocate space for sprintf() correctly.
#   0.92.24      - update MSTK to version 2.25 (updates to meshconvert, exoatt)
#   0.92.25      - update MSTK to version 2.26rc2 (adds fixes for pinchouts)
#
#   0.93.0       - defaulted to C++11, update Trilinos to 12.6.1
#   0.93.1       - update Boost to version 1.61.0

include(CMakeParseArguments)

MACRO(LIST_LENGTH var)
  SET(entries)
  FOREACH(e ${ARGN})
    SET(entries "${entries}.")
  ENDFOREACH(e)
  STRING(LENGTH "${entries}" ${var})
ENDMACRO(LIST_LENGTH)

# this macro appends version number defines to the tpl_versions.h include file
macro(amanzi_tpl_version_write)
  set(singleValueArgs FILENAME PREFIX)
  set(multiValueArgs VERSION)
  set(options "")
  
  cmake_parse_arguments(LOCAL "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  list_length(length ${LOCAL_VERSION})

  if (length GREATER 0) 
    list(GET LOCAL_VERSION 0 MAJOR)
    file(APPEND ${LOCAL_FILENAME} "#define ${LOCAL_PREFIX}_MAJOR ${MAJOR}\n")
  else()
    file(APPEND ${LOCAL_FILENAME} "#define ${LOCAL_PREFIX}_MAJOR\n")
  endif()

  if (length GREATER 1)
    list(GET LOCAL_VERSION 1 MINOR)
    file(APPEND ${LOCAL_FILENAME} "#define ${LOCAL_PREFIX}_MINOR ${MINOR}\n")
  else()
    file(APPEND ${LOCAL_FILENAME} "#define ${LOCAL_PREFIX}_MINOR\n")
  endif()

  if (length GREATER 2)
    list(GET LOCAL_VERSION 2 PATCH)
    file(APPEND ${LOCAL_FILENAME} "#define ${LOCAL_PREFIX}_PATCH ${PATCH}\n")
  else()
    file(APPEND ${LOCAL_FILENAME} "#define ${LOCAL_PREFIX}_PATCH\n")
  endif()

  file(APPEND ${LOCAL_FILENAME} "\n")

endmacro(amanzi_tpl_version_write)


set (AMANZI_TPLS_VERSION_MAJOR 0)
set (AMANZI_TPLS_VERSION_MINOR 93)
set (AMANZI_TPLS_VERSION_PATCH 1)
set (AMANZI_TPLS_VERSION ${AMANZI_TPLS_VERSION}.${AMANZI_TPLS_VERSION_MINOR}.${AMANZI_TPLS_VERSION_PATCH})
#   Not sure how to create a meaningful hash key for the collection

#
# TPL: Xerces
#
set(XERCES_VERSION_MAJOR 3)
set(XERCES_VERSION_MINOR 1)
set(XERCES_VERSION_PATCH 2)
set(XERCES_VERSION ${XERCES_VERSION_MAJOR}.${XERCES_VERSION_MINOR}.${XERCES_VERSION_PATCH})
set(XERCES_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(XERCES_ARCHIVE_FILE   xerces-c-${XERCES_VERSION}.tar.bz2)
set(XERCES_MD5_SUM        d987b8bb576aea456e92454781fe3615 ) 

#
# TPL: OpenMPI
#
set(OpenMPI_VERSION_MAJOR 1)
set(OpenMPI_VERSION_MINOR 4)
set(OpenMPI_VERSION_PATCH 4)
set(OpenMPI_VERSION ${OpenMPI_VERSION_MAJOR}.${OpenMPI_VERSION_MINOR}.${OpenMPI_VERSION_PATCH})
set(OpenMPI_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(OpenMPI_ARCHIVE_FILE   openmpi-${OpenMPI_VERSION}.tar.bz2)
set(OpenMPI_MD5_SUM        e58a1ea7b8af62453aaa0ddaee5f26a0) 

#
# TPL: CURL
#
set(CURL_VERSION_MAJOR 7)
set(CURL_VERSION_MINOR 37)
set(CURL_VERSION_PATCH 0)
set(CURL_VERSION ${CURL_VERSION_MAJOR}.${CURL_VERSION_MINOR}.${CURL_VERSION_PATCH})
set(CURL_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(CURL_ARCHIVE_FILE   curl-${CURL_VERSION}.tar.bz2)
set(CURL_MD5_SUM        7dda0cc2e4136f78d5801ac347be696b)

#
# TPL: zlib
#
set(ZLIB_VERSION_MAJOR 1)
set(ZLIB_VERSION_MINOR 2)
set(ZLIB_VERSION_PATCH 6)
set(ZLIB_VERSION ${ZLIB_VERSION_MAJOR}.${ZLIB_VERSION_MINOR}.${ZLIB_VERSION_PATCH})
set(ZLIB_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(ZLIB_ARCHIVE_FILE   zlib-${ZLIB_VERSION}.tar.gz)
set(ZLIB_MD5_SUM        618e944d7c7cd6521551e30b32322f4a) 

#
# TPL: METIS
#
set(METIS_VERSION_MAJOR 5)
set(METIS_VERSION_MINOR 1)
set(METIS_VERSION_PATCH 0)
set(METIS_VERSION ${METIS_VERSION_MAJOR}.${METIS_VERSION_MINOR}.${METIS_VERSION_PATCH})
set(METIS_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(METIS_ARCHIVE_FILE   metis-${METIS_VERSION}.tar.gz)
set(METIS_MD5_SUM        5465e67079419a69e0116de24fce58fe)

#
# TPL: CCSE
#
set(CCSE_VERSION_MAJOR 1)
set(CCSE_VERSION_MINOR 3)
set(CCSE_VERSION_PATCH 4)
set(CCSE_VERSION ${CCSE_VERSION_MAJOR}.${CCSE_VERSION_MINOR}.${CCSE_VERSION_PATCH})
set(AMANZI_DIR $ENV{AMANZI_DIR})
set(CCSE_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(CCSE_ARCHIVE_FILE   ccse-${CCSE_VERSION}.tar.gz) 
set(CCSE_MD5_SUM        faa52bb553cea8ca9ea436c1a7135b12)

#
# TPL: UnitTest
#
set(UnitTest_VERSION_MAJOR 1)
set(UnitTest_VERSION_MINOR 5)
set(UnitTest_VERSION ${UnitTest_VERSION_MAJOR}.${UnitTest_VERSION_MINOR})
set(UnitTest_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(UnitTest_ARCHIVE_FILE   unittest-cpp-${UnitTest_VERSION}.zip)
set(UnitTest_MD5_SUM      6f6e05fa07eeb2d44e5b11bd1f38865d) 

#
# TPL: Boost
#
set(Boost_VERSION_MAJOR 1)
set(Boost_VERSION_MINOR 61)
set(Boost_VERSION_PATCH 0)
set(Boost_VERSION        ${Boost_VERSION_MAJOR}.${Boost_VERSION_MINOR}.${Boost_VERSION_PATCH})
set(Boost_VERSION_STRING ${Boost_VERSION_MAJOR}_${Boost_VERSION_MINOR}_${Boost_VERSION_PATCH})
set(Boost_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(Boost_ARCHIVE_FILE   boost_${Boost_VERSION_STRING}.tar.bz2)
set(Boost_MD5_SUM        6095876341956f65f9d35939ccea1a9f)

#
# TPL: BoostCmake
#
set(BoostCmake_VERSION_MAJOR 1)
set(BoostCmake_VERSION_MINOR 46)
set(BoostCmake_VERSION_PATCH 1)
set(BoostCmake_VERSION        ${BoostCmake_VERSION_MAJOR}.${BoostCmake_VERSION_MINOR}.${BoostCmake_VERSION_PATCH})
set(BoostCmake_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(BoostCmake_ARCHIVE_FILE   boost-cmake-cmake-${BoostCmake_VERSION}.tar.gz)
set(BoostCmake_MD5_SUM        ) 

#
# TPL: HDF5
#
set(HDF5_VERSION_MAJOR 1)
set(HDF5_VERSION_MINOR 8)
set(HDF5_VERSION_PATCH 8)
set(HDF5_VERSION ${HDF5_VERSION_MAJOR}.${HDF5_VERSION_MINOR}.${HDF5_VERSION_PATCH})
set(HDF5_URL_STRING    "http://software.lanl.gov/ascem/tpls")
set(HDF5_ARCHIVE_FILE   hdf5-${HDF5_VERSION}.tar.gz)
set(HDF5_MD5_SUM        1196e668f5592bfb50d1de162eb16cff)      

#
# TPL: NetCDF
#
set(NetCDF_VERSION_MAJOR 4)
set(NetCDF_VERSION_MINOR 3)
set(NetCDF_VERSION_PATCH 2)
set(NetCDF_VERSION ${NetCDF_VERSION_MAJOR}.${NetCDF_VERSION_MINOR}.${NetCDF_VERSION_PATCH})
set(NetCDF_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(NetCDF_ARCHIVE_FILE   netcdf-${NetCDF_VERSION}.tar.gz)
set(NetCDF_MD5_SUM        2fd2365e1fe9685368cd6ab0ada532a0)

#
# TPL: NetCDF Fortran
#
set(NetCDF_Fortran_VERSION_MAJOR 4)
set(NetCDF_Fortran_VERSION_MINOR 2)
set(NetCDF_Fortran_VERSION ${NetCDF_Fortran_VERSION_MAJOR}.${NetCDF_Fortran_VERSION_MINOR})
set(NetCDF_Fortran_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(NetCDF_Fortran_ARCHIVE_FILE   netcdf-fortran-${NetCDF_Fortran_VERSION}.tar.gz)
set(NetCDF_Fortran_MD5_SUM        cc3bf530223e8f4aff93793b9f197bf3) 

#
# ASCEM-IO
#
set(ASCEMIO_VERSION_MAJOR 2)
set(ASCEMIO_VERSION_MINOR 2)
set(ASCEMIO_VERSION ${ASCEMIO_VERSION_MAJOR}.${ASCEMIO_VERSION_MINOR})
set(ASCEMIO_URL_STRING    "http://software.lanl.gov/ascem/tpls")
set(ASCEMIO_ARCHIVE_FILE   ascem-io-${ASCEMIO_VERSION}.tar.gz)
set(ASCEMIO_MD5_SUM       869820bacd4c289c8f320be58c1449a7)

#
# TPL: ExodusII
#
set(ExodusII_VERSION_MAJOR 6)
set(ExodusII_VERSION_MINOR 06)
set(ExodusII_VERSION ${ExodusII_VERSION_MAJOR}.${ExodusII_VERSION_MINOR})
set(ExodusII_URL_STRING    "http://software.lanl.gov/ascem/tpls")
set(ExodusII_ARCHIVE_FILE  exodus-${ExodusII_VERSION}.tar.gz)
set(ExodusII_MD5_SUM       cfd240dbc1251b08fb1d0ee2de40a44c)

#
# TPL: MSTK
#
set(MSTK_VERSION_MAJOR 2)
set(MSTK_VERSION_MINOR 26)
set(MSTK_VERSION_PATCH rc2)
set(MSTK_VERSION ${MSTK_VERSION_MAJOR}.${MSTK_VERSION_MINOR}${MSTK_VERSION_PATCH})
set(MSTK_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(MSTK_ARCHIVE_FILE   mstk-${MSTK_VERSION}.tgz)
set(MSTK_MD5_SUM        9063e949962c3ad6e16d1ce118e42bee)

#
# TPL: MOAB
#
set(MOAB_VERSION_MAJOR  r4276)
set(MOAB_VERSION_MINOR  )
set(MOAB_VERSION_PATCH  )
set(MOAB_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(MOAB_ARCHIVE_FILE   MOAB-${MOAB_VERSION}.tar.gz)
set(MOAB_MD5_SUM        49da04e8905f6d730d92521e7ca7400e) 

#
# TPL: HYPRE
#
set(HYPRE_VERSION_MAJOR  2)
set(HYPRE_VERSION_MINOR  10)
set(HYPRE_VERSION_PATCH  0b)
set(HYPRE_VERSION  ${HYPRE_VERSION_MAJOR}.${HYPRE_VERSION_MINOR}.${HYPRE_VERSION_PATCH})
set(HYPRE_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(HYPRE_ARCHIVE_FILE   hypre-${HYPRE_VERSION}.tar.gz)
set(HYPRE_MD5_SUM        768be38793a35bb5d055905b271f5b8e) 

#
# TPL: ParMetis
#
set(ParMetis_VERSION_MAJOR  4)
set(ParMetis_VERSION_MINOR  0)
set(ParMetis_VERSION_PATCH  3a)
set(ParMetis_VERSION  ${ParMetis_VERSION_MAJOR}.${ParMetis_VERSION_MINOR}.${ParMetis_VERSION_PATCH})
set(ParMetis_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(ParMetis_ARCHIVE_FILE   parmetis-${ParMetis_VERSION}.tar.gz)
set(ParMetis_MD5_SUM        56ac6ebf6e7e8a522fa053c799dc7a92)

#
# TPL: SuperLU (Built by PETSc!)
#
set(SuperLU_VERSION_MAJOR  4)
set(SuperLU_VERSION_MINOR  3)
set(SuperLU_VERSION  ${SuperLU_VERSION_MAJOR}.${SuperLU_VERSION_MINOR})
set(SuperLU_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(SuperLU_ARCHIVE_FILE   superlu_${SuperLU_VERSION}.tar.gz)
set(SuperLU_MD5_SUM        b72c6309f25e9660133007b82621ba7c)

#
# TPL: SuperLU Distrib (Built by PETSc!)
#
set(SuperLUDist_VERSION_MAJOR  3)
set(SuperLUDist_VERSION_MINOR  3)
set(SuperLUDist_VERSION  ${SuperLUDist_VERSION_MAJOR}.${SuperLUDist_VERSION_MINOR})
set(SuperLUDist_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(SuperLUDist_ARCHIVE_FILE   superlu_dist_${SuperLUDist_VERSION}.tar.gz)
set(SuperLUDist_MD5_SUM        b72c6309f25e9660133007b82621ba7c)

#
# TPL: PETSc
#
set(PETSc_VERSION_MAJOR  3)
set(PETSc_VERSION_MINOR  5)
set(PETSc_VERSION_PATCH  2)
set(PETSc_VERSION  ${PETSc_VERSION_MAJOR}.${PETSc_VERSION_MINOR}.${PETSc_VERSION_PATCH})
set(PETSc_ARCHIVE_VERSION ${PETSc_VERSION_MAJOR}.${PETSc_VERSION_MINOR}.${PETSc_VERSION_PATCH})
set(PETSc_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(PETSc_ARCHIVE_FILE   petsc-${PETSc_ARCHIVE_VERSION}.tar.gz)
set(PETSc_MD5_SUM        ad170802b3b058b5deb9cd1f968e7e13)

#
# TPL: Trilinos
#
set(Trilinos_VERSION_MAJOR 12)
set(Trilinos_VERSION_MINOR 6)
set(Trilinos_VERSION_PATCH 3)
set(Trilinos_VERSION ${Trilinos_VERSION_MAJOR}-${Trilinos_VERSION_MINOR}-${Trilinos_VERSION_PATCH})
set(Trilinos_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(Trilinos_ARCHIVE_FILE   Trilinos-trilinos-release-${Trilinos_VERSION}.tar.gz)
set(Trilinos_MD5_SUM        8de5cc00981a0ca0defea6199b2fe4c1)

#
# TPL: SEACAS
#  SEACAS is available in Trilinos 10.8 and above
set(SEACAS_VERSION_MAJOR 12)
set(SEACAS_VERSION_MINOR 6)
set(SEACAS_VERSION_PATCH 3)
set(SEACAS_VERSION ${SEACAS_VERSION_MAJOR}-${SEACAS_VERSION_MINOR}-${SEACAS_VERSION_PATCH})
set(SEACAS_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(SEACAS_ARCHIVE_FILE   Trilinos-trilinos-release-${SEACAS_VERSION}.tar.gz)
set(SEACAS_MD5_SUM        8de5cc00981a0ca0defea6199b2fe4c1)

#
# TPL: PFlotran
#
set(PFLOTRAN_VERSION_MAJOR 0)
set(PFLOTRAN_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(PFLOTRAN_ARCHIVE_FILE   pflotran-dev-611092f80ddb.tar.gz)
set(PFLOTRAN_MD5_SUM        e18997dd7de5523c9bef8489a0a2dd24)

#
# TPL: Alquimia
#
set(ALQUIMIA_VERSION_MAJOR 1)
set(ALQUIMIA_VERSION_MINOR 0)
set(ALQUIMIA_VERSION_PATCH 3)
set(ALQUIMIA_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(ALQUIMIA_ARCHIVE_FILE   alquimia-dev-1.0.3.tar.gz)
set(ALQUIMIA_MD5_SUM        13525a9d44df905fed45da75eca6fedb)


#
# TPL: Silo
#
set(Silo_VERSION_MAJOR 4)
set(Silo_VERSION_MINOR 10)
set(Silo_VERSION_PATCH 2)
set(Silo_VERSION  ${Silo_VERSION_MAJOR}.${Silo_VERSION_MINOR}.${Silo_VERSION_PATCH})
set(Silo_URL_STRING "http://software.lanl.gov/ascem/tpls")
set(Silo_ARCHIVE_FILE silo-4.10.2.tar.gz)
set(Silo_MD5_SUM 9ceac777a2f2469ac8cef40f4fab49c8)

