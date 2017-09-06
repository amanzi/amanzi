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
#   0.93.2       - update Alquimia to version 1.0.3
#   0.93.3       - update Alquimia to version 1.0.4

#   0.94.1       - updates several TPLs, new versions are:
#                - Trilinos 12.10.1
#                - zlib 1.2.11
#                - hdf5 1.8.12
#                - netcdf 4.4.1.1
#                - netcdf-fortran 4.4.4
#                - boost 1.6.3
#   0.94.2       - update MTSK to version 2.28rc1
#   0.94.3       - update MSTK to version 3.00 (incompatible - need to update #defines)
#   0.94.4       - update MSTK to version 3.01
#   0.94.5       - restored Alquimia to version 1.0.4
#   0.94.6       - Added CrunchTope package, hash version c31ecb9
#   0.94.7       - update UnitTest++ to version 2.0.0
#                - update Hypre to version 2.11.2    
#   0.94.8       - removed ExodusII as independent TPL  
#   0.94.9       - update PFloTran to version dev-c8df814cb6fa
#                - update PETSc to xsdk-0.2.0 (native 3.7.5)
#                - update SuperLU to 5.2.1
#                - update SuperLU_dist to xsdk-0.2.0 (native 5.1.3)
#                - update Alquimia to xsdk-0.2.0 (native 1.0.4)
#                - update Hypre to xsdk-0.2.0 (native 2.11.2)

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


#
# TPLs and XSDK versions 
#
set(AMANZI_TPLS_VERSION_MAJOR 0)
set(AMANZI_TPLS_VERSION_MINOR 94)
set(AMANZI_TPLS_VERSION_PATCH 9)
set(AMANZI_TPLS_VERSION ${AMANZI_TPLS_VERSION_MAJOR}.${AMANZI_TPLS_VERSION_MINOR}.${AMANZI_TPLS_VERSION_PATCH})
# Not sure how to create a meaningful hash key for the collection

set(XSDK_VERSION "0.2.0")

#
# Default location on GitHub
#
set (AMANZI_TPLS_DOWNLOAD_URL "https://raw.githubusercontent.com/amanzi/amanzi-tpls/master/src")

#
# TPL: Xerces
#
set(XERCES_VERSION_MAJOR 3)
set(XERCES_VERSION_MINOR 1)
set(XERCES_VERSION_PATCH 2)
set(XERCES_VERSION ${XERCES_VERSION_MAJOR}.${XERCES_VERSION_MINOR}.${XERCES_VERSION_PATCH})
set(XERCES_URL_STRING     ${AMANZI_TPLS_DOWNLOAD_URL})
set(XERCES_ARCHIVE_FILE   xerces-c-${XERCES_VERSION}.tar.bz2)
set(XERCES_MD5_SUM        d987b8bb576aea456e92454781fe3615 ) 

#
# TPL: OpenMPI
#
set(OpenMPI_VERSION_MAJOR 1)
set(OpenMPI_VERSION_MINOR 4)
set(OpenMPI_VERSION_PATCH 4)
set(OpenMPI_VERSION ${OpenMPI_VERSION_MAJOR}.${OpenMPI_VERSION_MINOR}.${OpenMPI_VERSION_PATCH})
set(OpenMPI_URL_STRING     ${ASCEM_TPLS_DOWNLOAD_URL})
set(OpenMPI_ARCHIVE_FILE   openmpi-${OpenMPI_VERSION}.tar.bz2)
set(OpenMPI_MD5_SUM        e58a1ea7b8af62453aaa0ddaee5f26a0) 

#
# TPL: CURL
#
set(CURL_VERSION_MAJOR 7)
set(CURL_VERSION_MINOR 37)
set(CURL_VERSION_PATCH 0)
set(CURL_VERSION ${CURL_VERSION_MAJOR}.${CURL_VERSION_MINOR}.${CURL_VERSION_PATCH})
set(CURL_URL_STRING     ${AMANZI_TPLS_DOWNLOAD_URL})
set(CURL_ARCHIVE_FILE   curl-${CURL_VERSION}.tar.bz2)
set(CURL_MD5_SUM        7dda0cc2e4136f78d5801ac347be696b)

#
# TPL: zlib
#
set(ZLIB_VERSION_MAJOR 1)
set(ZLIB_VERSION_MINOR 2)
set(ZLIB_VERSION_PATCH 11)
set(ZLIB_VERSION ${ZLIB_VERSION_MAJOR}.${ZLIB_VERSION_MINOR}.${ZLIB_VERSION_PATCH})
set(ZLIB_URL_STRING     ${AMANZI_TPLS_DOWNLOAD_URL})
set(ZLIB_ARCHIVE_FILE   zlib-${ZLIB_VERSION}.tar.gz)
set(ZLIB_MD5_SUM        1c9f62f0778697a09d36121ead88e08e) 

#
# TPL: METIS
#
set(METIS_VERSION_MAJOR 5)
set(METIS_VERSION_MINOR 1)
set(METIS_VERSION_PATCH 0)
set(METIS_VERSION ${METIS_VERSION_MAJOR}.${METIS_VERSION_MINOR}.${METIS_VERSION_PATCH})
set(METIS_URL_STRING     ${AMANZI_TPLS_DOWNLOAD_URL})
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
set(CCSE_URL_STRING     ${AMANZI_TPLS_DOWNLOAD_URL})
set(CCSE_ARCHIVE_FILE   ccse-${CCSE_VERSION}.tar.gz) 
set(CCSE_MD5_SUM        faa52bb553cea8ca9ea436c1a7135b12)

#
# TPL: UnitTest
#
set(UnitTest_VERSION_MAJOR 2)
set(UnitTest_VERSION_MINOR 0)
set(UnitTest_VERSION_PATCH 0)
set(UnitTest_VERSION ${UnitTest_VERSION_MAJOR}.${UnitTest_VERSION_MINOR}.${UnitTest_VERSION_PATCH})
set(UnitTest_URL_STRING     ${AMANZI_TPLS_DOWNLOAD_URL})
set(UnitTest_ARCHIVE_FILE   unittest-cpp-${UnitTest_VERSION}.tgz)
set(UnitTest_MD5_SUM      29f958e355e516e7ab016b467974728d) 

#
# TPL: Boost
#
set(Boost_VERSION_MAJOR 1)
set(Boost_VERSION_MINOR 63)
set(Boost_VERSION_PATCH 0)
set(Boost_VERSION        ${Boost_VERSION_MAJOR}.${Boost_VERSION_MINOR}.${Boost_VERSION_PATCH})
set(Boost_VERSION_STRING ${Boost_VERSION_MAJOR}_${Boost_VERSION_MINOR}_${Boost_VERSION_PATCH})
set(Boost_URL_STRING     ${AMANZI_TPLS_DOWNLOAD_URL})
set(Boost_ARCHIVE_FILE   boost_${Boost_VERSION_STRING}.tar.bz2)
set(Boost_MD5_SUM        1c837ecd990bb022d07e7aab32b09847)

#
# TPL: HDF5
#
set(HDF5_VERSION_MAJOR 1)
set(HDF5_VERSION_MINOR 8)
set(HDF5_VERSION_PATCH 18)
set(HDF5_VERSION ${HDF5_VERSION_MAJOR}.${HDF5_VERSION_MINOR}.${HDF5_VERSION_PATCH})
set(HDF5_URL_STRING     ${AMANZI_TPLS_DOWNLOAD_URL})
set(HDF5_ARCHIVE_FILE   hdf5-${HDF5_VERSION}.tar.gz)
set(HDF5_MD5_SUM        dd2148b740713ca0295442ec683d7b1c)


#
# TPL: NetCDF
#
set(NetCDF_VERSION_MAJOR 4)
set(NetCDF_VERSION_MINOR 4)
set(NetCDF_VERSION_PATCH 1.1)
set(NetCDF_VERSION ${NetCDF_VERSION_MAJOR}.${NetCDF_VERSION_MINOR}.${NetCDF_VERSION_PATCH})
set(NetCDF_URL_STRING     ${AMANZI_TPLS_DOWNLOAD_URL})
set(NetCDF_ARCHIVE_FILE   netcdf-${NetCDF_VERSION}.tar.gz)
set(NetCDF_MD5_SUM        503a2d6b6035d116ed53b1d80c811bda)

#
# TPL: NetCDF Fortran
#
set(NetCDF_Fortran_VERSION_MAJOR 4)
set(NetCDF_Fortran_VERSION_MINOR 2)
set(NetCDF_Fortran_VERSION ${NetCDF_Fortran_VERSION_MAJOR}.${NetCDF_Fortran_VERSION_MINOR})
set(NetCDF_Fortran_URL_STRING     ${AMANZI_TPLS_DOWNLOAD_URL})
set(NetCDF_Fortran_ARCHIVE_FILE   netcdf-fortran-${NetCDF_Fortran_VERSION}.tar.gz)
set(NetCDF_Fortran_MD5_SUM        cc3bf530223e8f4aff93793b9f197bf3) 

#
# ASCEM-IO
#
set(ASCEMIO_VERSION_MAJOR 2)
set(ASCEMIO_VERSION_MINOR 2)
set(ASCEMIO_VERSION ${ASCEMIO_VERSION_MAJOR}.${ASCEMIO_VERSION_MINOR})
set(ASCEMIO_URL_STRING    ${AMANZI_TPLS_DOWNLOAD_URL})
set(ASCEMIO_ARCHIVE_FILE   ascem-io-${ASCEMIO_VERSION}.tar.gz)
set(ASCEMIO_MD5_SUM       869820bacd4c289c8f320be58c1449a7)

#
# TPL: MSTK
#
set(MSTK_VERSION_MAJOR 3)
set(MSTK_VERSION_MINOR 0)
set(MSTK_VERSION_PATCH 1)
set(MSTK_VERSION ${MSTK_VERSION_MAJOR}.${MSTK_VERSION_MINOR}.${MSTK_VERSION_PATCH})
set(MSTK_URL_STRING     ${AMANZI_TPLS_DOWNLOAD_URL})
set(MSTK_ARCHIVE_FILE   mstk-${MSTK_VERSION}.tar.gz)
set(MSTK_MD5_SUM        d0761665844b1f956ef7cb3d80d68d88)

#
# TPL: MOAB
#
set(MOAB_VERSION_MAJOR  5)
set(MOAB_VERSION_MINOR  0)
set(MOAB_VERSION_PATCH  0)
set(MOAB_VERSION ${MOAB_VERSION_MAJOR}.${MOAB_VERSION_MINOR}.${MOAB_VERSION_PATCH})
set(MOAB_URL_STRING     ftp://ftp.mcs.anl.gov/pub/fathom)
set(MOAB_ARCHIVE_FILE   moab-${MOAB_VERSION}.tar.gz)
set(MOAB_MD5_SUM        1840ca02366f4d3237d44af63e239e3b) 

#
# TPL: HYPRE
#
set(HYPRE_VERSION_MAJOR  2)
set(HYPRE_VERSION_MINOR  11)
set(HYPRE_VERSION_PATCH  2)
set(HYPRE_VERSION  ${HYPRE_VERSION_MAJOR}.${HYPRE_VERSION_MINOR}.${HYPRE_VERSION_PATCH})
set(HYPRE_URL_STRING     "https://github.com/LLNL/hypre/archive/")
set(HYPRE_ARCHIVE_FILE   xsdk-${XSDK_VERSION}.tar.gz)
set(HYPRE_SAVEAS_FILE    hypre-${HYPRE_VERSION}.tar.gz)
set(HYPRE_MD5_SUM        fc9474058560602e9be2ce618db7fd14) 

#
# TPL: ParMetis
#
set(ParMetis_VERSION_MAJOR  4)
set(ParMetis_VERSION_MINOR  0)
set(ParMetis_VERSION_PATCH  3a)
set(ParMetis_VERSION  ${ParMetis_VERSION_MAJOR}.${ParMetis_VERSION_MINOR}.${ParMetis_VERSION_PATCH})
set(ParMetis_URL_STRING     ${AMANZI_TPLS_DOWNLOAD_URL})
set(ParMetis_ARCHIVE_FILE   parmetis-${ParMetis_VERSION}.tar.gz)
set(ParMetis_MD5_SUM        56ac6ebf6e7e8a522fa053c799dc7a92)

#
# TPL: SuperLU (Built by PETSc!)
#
set(SuperLU_VERSION_MAJOR  5)
set(SuperLU_VERSION_MINOR  2)
set(SuperLU_VERSION_PATCH  1)
set(SuperLU_VERSION  ${SuperLU_VERSION_MAJOR}.${SuperLU_VERSION_MINOR}.${SuperLU_VERSION_PATCH})
set(SuperLU_URL_STRING     "http://crd-legacy.lbl.gov/~xiaoye/SuperLU")
set(SuperLU_ARCHIVE_FILE   superlu_${SuperLU_VERSION}.tar.gz)
set(SuperLU_SAVEAS_FILE    superlu_${SuperLU_VERSION}.tar.gz)
set(SuperLU_MD5_SUM        3a1a9bff20cb06b7d97c46d337504447)

#
# TPL: SuperLU Distrib (Built by PETSc!)
#
set(SuperLUDist_VERSION_MAJOR  5)
set(SuperLUDist_VERSION_MINOR  1)
set(SuperLUDist_VERSION_PATCH  3)
set(SuperLUDist_VERSION  ${SuperLUDist_VERSION_MAJOR}.${SuperLUDist_VERSION_MINOR}.${SuperLUDist_VERSION_PATCH})
set(SuperLUDist_URL_STRING     "https://github.com/xiaoyeli/superlu_dist/archive")
set(SuperLUDist_ARCHIVE_FILE   xsdk-${XSDK_VERSION}.tar.gz)
set(SuperLUDist_SAVEAS_FILE    superlu_dist_${SuperLUDist_VERSION}.tar.gz)
set(SuperLUDist_MD5_SUM        9ccd1915dd06f167ed8dca7b14bbcedb)

#
# TPL: PETSc
#
set(PETSc_VERSION_MAJOR  3)
set(PETSc_VERSION_MINOR  7)
set(PETSc_VERSION_PATCH  5)
set(PETSc_VERSION  ${PETSc_VERSION_MAJOR}.${PETSc_VERSION_MINOR}.${PETSc_VERSION_PATCH})
set(PETSc_ARCHIVE_VERSION ${PETSc_VERSION_MAJOR}.${PETSc_VERSION_MINOR}.${PETSc_VERSION_PATCH})
set(PETSc_URL_STRING     "https://bitbucket.org/petsc/petsc/get")
set(PETSc_ARCHIVE_FILE   xsdk-${XSDK_VERSION}.tar.gz)
set(PETSc_SAVEAS_FILE    petsc-${PETSc_ARCHIVE_VERSION}.tar.gz)
set(PETSc_MD5_SUM        41a10be8bbf9d13f137873a2d52c6715)

#
# TPL: Trilinos
#
set(Trilinos_VERSION_MAJOR 12)
set(Trilinos_VERSION_MINOR 10)
set(Trilinos_VERSION_PATCH 1)
set(Trilinos_VERSION ${Trilinos_VERSION_MAJOR}-${Trilinos_VERSION_MINOR}-${Trilinos_VERSION_PATCH})
set(Trilinos_URL_STRING     "https://github.com/trilinos/Trilinos/archive")
set(Trilinos_ARCHIVE_FILE   trilinos-release-${Trilinos_VERSION}.tar.gz)
set(Trilinos_MD5_SUM        40f28628b63310f9bd17c26d9ebe32b1)

#
# TPL: SEACAS
#
set(SEACAS_VERSION_MAJOR 173a1e6)
set(SEACAS_VERSION_MINOR 0)
set(SEACAS_VERSION_PATCH 0)
set(SEACAS_VERSION ${SEACAS_VERSION_MAJOR})
set(SEACAS_URL_STRING     ${AMANZI_TPLS_DOWNLOAD_URL})
set(SEACAS_ARCHIVE_FILE   seacas-${SEACAS_VERSION}.tgz)
set(SEACAS_MD5_SUM        3235d1b885ee8e1a04408382f50bd0f0)

#
# TPL: PFlotran
#
set(PFLOTRAN_VERSION_MAJOR 0)
set(PFLOTRAN_VERSION_MINOR 2)
set(PFLOTRAN_VERSION_PATCH 0)
set(PFLOTRAN_VERSION ${PFLOTRAN_VERSION_MAJOR}.${PFLOTRAN_VERSION_MINOR}.${PFLOTRAN_VERSION_PATCH})
set(PFLOTRAN_URL_STRING     "https://bitbucket.org/pflotran/pflotran/get")
set(PFLOTRAN_ARCHIVE_FILE   xsdk-${XSDK_VERSION}-rc2.tar.gz)
set(PFLOTRAN_SAVEAS_FILE    pflotran-${PFLOTRAN_VERSION}.tar.gz)
set(PFLOTRAN_MD5_SUM        80a214c394bbd4230c2ddc0ba177c8ea)

#
# TPL: Alquimia
#
set(ALQUIMIA_VERSION_MAJOR 1)
set(ALQUIMIA_VERSION_MINOR 0)
set(ALQUIMIA_VERSION_PATCH 4)
set(ALQUIMIA_VERSION ${ALQUIMIA_VERSION_MAJOR}.${ALQUIMIA_VERSION_MINOR}.${ALQUIMIA_VERSION_PATCH})
set(ALQUIMIA_URL_STRING     https://github.com/LBL-EESA/alquimia-dev/archive)
set(ALQUIMIA_ARCHIVE_FILE   xsdk-${XSDK_VERSION}.tar.gz)
set(ALQUIMIA_SAVEAS_FILE    alquimia-dev-${ALQUIMIA_VERSION}.tar.gz)
set(ALQUIMIA_MD5_SUM        c9ad5100a4f064c3cdc49bcebd0e508e)

#
# TPL: Silo
#
set(Silo_VERSION_MAJOR 4)
set(Silo_VERSION_MINOR 10)
set(Silo_VERSION_PATCH 2)
set(Silo_VERSION  ${Silo_VERSION_MAJOR}.${Silo_VERSION_MINOR}.${Silo_VERSION_PATCH})
set(Silo_URL_STRING "https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2")
set(Silo_ARCHIVE_FILE silo-4.10.2.tar.gz)
set(Silo_MD5_SUM 9ceac777a2f2469ac8cef40f4fab49c8)

#
# TPL: CrunchTope
#
set(CRUNCHTOPE_VERSION_MAJOR 160915)
set(CRUNCHTOPE_VERSION_MINOR c31ecb9)
set(CRUNCHTOPE_VERSION_PATCH 0)
set(CRUNCHTOPE_VERSION  ${CRUNCHTOPE_VERSION_MAJOR}.${CRUNCHTOPE_VERSION_MINOR}.${CRUNCHTOPE_VERSION_PATCH})
set(CRUNCHTOPE_URL_STRING ${AMANZI_TPLS_DOWNLOAD_URL})
set(CRUNCHTOPE_ARCHIVE_FILE CrunchTope_160915-c31ecb9.tgz)
set(CRUNCHTOPE_MD5_SUM 84c38ca70da8f0e14cce3841dbbb4c0b)

