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
#   0.91.10      - updated MSTK to verion 2.10rc5
#

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
set (AMANZI_TPLS_VERSION_MINOR 91)
set (AMANZI_TPLS_VERSION_PATCH 10)
set (AMANZI_TPLS_VERSION ${AMANZI_TPLS_VERSION}.${AMANZI_TPLS_VERSION_MINOR}.${AMANZI_TPLS_VERSION_PATCH})
#   Not sure how to create a meaningful hash key for the collection

#
# TPL: Xerces
#
set(XERCES_VERSION_MAJOR 3)
set(XERCES_VERSION_MINOR 1)
set(XERCES_VERSION_PATCH 1)
set(XERCES_VERSION ${XERCES_VERSION_MAJOR}.${XERCES_VERSION_MINOR}.${XERCES_VERSION_PATCH})
set(XERCES_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(XERCES_ARCHIVE_FILE   xerces-c-${XERCES_VERSION}.tar.gz)
set(XERCES_MD5_SUM        6a8ec45d83c8cfb1584c5a5345cb51ae ) 

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
set(CURL_VERSION_MINOR 21)
set(CURL_VERSION_PATCH 6)
set(CURL_VERSION ${CURL_VERSION_MAJOR}.${CURL_VERSION_MINOR}.${CURL_VERSION_PATCH})
set(CURL_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(CURL_ARCHIVE_FILE   curl-${CURL_VERSION}.tar.gz)
set(CURL_MD5_SUM        c502b67898b4a1bd687fe1b86419a44b) 

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
set(METIS_VERSION_MAJOR 4)
set(METIS_VERSION_MINOR 0)
set(METIS_VERSION_PATCH 3)
set(METIS_VERSION ${METIS_VERSION_MAJOR}.${METIS_VERSION_MINOR}.${METIS_VERSION_PATCH})
set(METIS_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(METIS_ARCHIVE_FILE   metis-${METIS_VERSION}.tar.gz)
set(METIS_MD5_SUM        d3848b454532ef18dc83e4fb160d1e10) 

#
# TPL: CCSE
#
set(CCSE_VERSION_MAJOR 1)
set(CCSE_VERSION_MINOR 2)
set(CCSE_VERSION_PATCH 4)
set(CCSE_VERSION ${CCSE_VERSION_MAJOR}.${CCSE_VERSION_MINOR}.${CCSE_VERSION_PATCH})
set(AMANZI_DIR $ENV{AMANZI_DIR})
set(CCSE_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(CCSE_ARCHIVE_FILE   ccse-${CCSE_VERSION}.tar.gz) 
set(CCSE_MD5_SUM        4a565639347a9bab5f6db70d89dfbb8b)

#
# TPL: UnitTest
#
set(UnitTest_VERSION_MAJOR 1)
set(UnitTest_VERSION_MINOR 4)
set(UnitTest_VERSION ${UnitTest_VERSION_MAJOR}.${UnitTest_VERSION_MINOR})
set(UnitTest_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(UnitTest_ARCHIVE_FILE   unittest-cpp-${UnitTest_VERSION}.zip)
set(UnitTest_MD5_SUM       bd373a53403ed51ea1bbb60b1952d7e3) 

#
# TPL: Boost
#
set(Boost_VERSION_MAJOR 1)
set(Boost_VERSION_MINOR 51)
set(Boost_VERSION_PATCH 0)
set(Boost_VERSION        ${Boost_VERSION_MAJOR}.${Boost_VERSION_MINOR}.${Boost_VERSION_PATCH})
set(Boost_VERSION_STRING ${Boost_VERSION_MAJOR}_${Boost_VERSION_MINOR}_${Boost_VERSION_PATCH})
set(Boost_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(Boost_ARCHIVE_FILE   boost_${Boost_VERSION_STRING}.tar.bz2)
set(Boost_MD5_SUM        4b6bd483b692fd138aef84ed2c8eb679)

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
set(NetCDF_VERSION_MINOR 2)
set(NetCDF_VERSION_PATCH 1.1)
set(NetCDF_VERSION ${NetCDF_VERSION_MAJOR}.${NetCDF_VERSION_MINOR}.${NetCDF_VERSION_PATCH})
set(NetCDF_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(NetCDF_ARCHIVE_FILE   netcdf-${NetCDF_VERSION}.tar.gz)
set(NetCDF_MD5_SUM        5eebcf19e6ac78a61c73464713cbfafc)

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
option(ENABLE_EXODUS498 "Use the 4.98 version of ExodusII" OFF)
if ( ENABLE_EXODUS498 )
  set(ExodusII_VERSION_MAJOR 4)
  set(ExodusII_VERSION_MINOR 98)
  set(ExodusII_MD5_SUM        4480e641d6ada58f5d8ecb7172e76791) 
else()
  set(ExodusII_VERSION_MAJOR 5)
  set(ExodusII_VERSION_MINOR 22)
  set(ExodusII_MD5_SUM        b5f537d7028f2c8cf6b4ba3e8a469dbe) 
endif()
set(ExodusII_VERSION ${ExodusII_VERSION_MAJOR}.${ExodusII_VERSION_MINOR})
set(ExodusII_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(ExodusII_ARCHIVE_FILE   exodusii-${ExodusII_VERSION}.tar.gz)

#
# TPL: MSTK
#
set(MSTK_VERSION_MAJOR 2)
set(MSTK_VERSION_MINOR 10)
set(MSTK_VERSION_PATCH rc5)
set(MSTK_VERSION ${MSTK_VERSION_MAJOR}.${MSTK_VERSION_MINOR}${MSTK_VERSION_PATCH})
set(MSTK_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(MSTK_ARCHIVE_FILE   mstk-${MSTK_VERSION}.tgz)
set(MSTK_MD5_SUM        9baf10c1cb3797d23d1d460bd39d66a8)    

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
set(HYPRE_VERSION_MINOR  8)
set(HYPRE_VERSION_PATCH  0b)
set(HYPRE_VERSION  ${HYPRE_VERSION_MAJOR}.${HYPRE_VERSION_MINOR}.${HYPRE_VERSION_PATCH})
set(HYPRE_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(HYPRE_ARCHIVE_FILE   hypre-${HYPRE_VERSION}.tar.gz)
set(HYPRE_MD5_SUM        6b4db576c68d2072e48efbc00ea58489) 

#
# TPL: ParMetis (Built by PETSc!)
#
set(ParMetis_VERSION_MAJOR  3)
set(ParMetis_VERSION_MINOR  2)
set(ParMetis_VERSION_PATCH  0-p1)
set(ParMetis_VERSION  ${ParMetis_VERSION_MAJOR}.${ParMetis_VERSION_MINOR}.${ParMetis_VERSION_PATCH})
set(ParMetis_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(ParMetis_ARCHIVE_FILE   ParMetis-${ParMetis_VERSION}.tar.gz)
set(ParMetis_MD5_SUM        f17ec2aeacc04f67f8b69f28cae4079f) 

#
# TPL: SuperLU (Built by PETSc!)
#
set(SuperLU_VERSION_MAJOR  4)
set(SuperLU_VERSION_MINOR  2)
set(SuperLU_VERSION  ${SuperLU_VERSION_MAJOR}.${SuperLU_VERSION_MINOR})
set(SuperLU_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(SuperLU_ARCHIVE_FILE   superlu_${SuperLU_VERSION}.tar.gz)
set(SuperLU_MD5_SUM        565602cf69e425874c2525f8b96e9bb1)

#
# TPL: SuperLU Distrib (Built by PETSc!)
#
set(SuperLUDist_VERSION_MAJOR  2)
set(SuperLUDist_VERSION_MINOR  5)
set(SuperLUDist_VERSION  ${SuperLUDist_VERSION_MAJOR}.${SuperLUDist_VERSION_MINOR})
set(SuperLUDist_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(SuperLUDist_ARCHIVE_FILE   superlu_dist_${SuperLUDist_VERSION}.tar.gz)
set(SuperLUDist_MD5_SUM        2194ae8f9786e396a721cf4d41045566)
 
#
# TPL: PETSc
#
set(PETSc_VERSION_MAJOR  3)
set(PETSc_VERSION_MINOR  2)
set(PETSc_VERSION_PATCH  7)
set(PETSc_VERSION  ${PETSc_VERSION_MAJOR}.${PETSc_VERSION_MINOR}.${PETSc_VERSION_PATCH})
set(PETSc_ARCHIVE_VERSION ${PETSc_VERSION_MAJOR}.${PETSc_VERSION_MINOR}-p${PETSc_VERSION_PATCH})
set(PETSc_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(PETSc_ARCHIVE_FILE   petsc-${PETSc_ARCHIVE_VERSION}.tar.gz)
set(PETSc_MD5_SUM        b9b5b42ffb6c619e4f7ee6b29134dc5f)

#
# TPL: Trilinos
#
set(Trilinos_VERSION_MAJOR 11)
set(Trilinos_VERSION_MINOR 0)
set(Trilinos_VERSION_PATCH 3)
set(Trilinos_VERSION ${Trilinos_VERSION_MAJOR}.${Trilinos_VERSION_MINOR}.${Trilinos_VERSION_PATCH})
set(Trilinos_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(Trilinos_ARCHIVE_FILE   trilinos-${Trilinos_VERSION}-Source.tar.bz2)
set(Trilinos_MD5_SUM        edd4cbb84ff5874dc96e846ceeab50b7) 

#
# TPL: SEACAS
#  SEACAS is available in Trilinos 10.8 and above
set(SEACAS_VERSION_MAJOR 11)
set(SEACAS_VERSION_MINOR 0)
set(SEACAS_VERSION_PATCH 3)
set(SEACAS_VERSION ${SEACAS_VERSION_MAJOR}.${SEACAS_VERSION_MINOR}.${SEACAS_VERSION_PATCH})
set(SEACAS_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(SEACAS_ARCHIVE_FILE   trilinos-${SEACAS_VERSION}-Source.tar.bz2)
set(SEACAS_MD5_SUM        edd4cbb84ff5874dc96e846ceeab50b7) 

#
# TPL: PFlotran
#
set(PFLOTRAN_VERSION_MAJOR 0)
set(PFLOTRAN_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(PFLOTRAN_ARCHIVE_FILE   pflotran-dev-07032013-clean.tgz)
set(PFLOTRAN_MD5_SUM        f81f285837e4490992e9dda2159e1baa)

#
# TPL: Alquimia
#
set(ALQUIMIA_VERSION_MAJOR 0)
set(ALQUIMIA_VERSION_MINOR 1)
set(ALQUIMIA_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(ALQUIMIA_ARCHIVE_FILE   alquimia-0.1.tar.gz)
set(ALQUIMIA_MD5_SUM        8001a8721a002937801562bcc242824f)
