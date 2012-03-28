#  -*- mode: cmake -*-

#
# TPLVersions
#    Define the versions, approved download locations for each TPL
#
#

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
set(ZLIB_URL_STRING     "http://zlib.net")
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
set(CCSE_VERSION_MAJOR 0)
set(CCSE_VERSION_MINOR 1)
set(CCSE_VERSION_PATCH 10)
set(CCSE_VERSION ${CCSE_VERSION_MAJOR}.${CCSE_VERSION_MINOR}.${CCSE_VERSION_PATCH})
set(CCSE_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(CCSE_ARCHIVE_FILE   ccse-${CCSE_VERSION}.tar.gz)
set(CCSE_MD5_SUM        0914c9ef955e64490335a54d583dd648) 


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
# TPL: CGNS
#
set(CGNS_VERSION_MAJOR 2)
set(CGNS_VERSION_MINOR 5)
set(CGNS_VERSION_PATCH 4)
set(CGNS_VERSION ${CGNS_VERSION_MAJOR}.${CGNS_VERSION_MINOR}-${CGNS_VERSION_PATCH})
set(CGNS_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(CGNS_ARCHIVE_FILE   cgnslib_${CGNS_VERSION}.tar.gz)
set(CGNS_MD5_SUM        42063efdf726c81300a51c3495d3224e) 


#
# TPL: Boost
#
set(Boost_VERSION_MAJOR 1)
set(Boost_VERSION_MINOR 48)
set(Boost_VERSION_PATCH 0)
set(Boost_VERSION        ${Boost_VERSION_MAJOR}.${Boost_VERSION_MINOR}.${Boost_VERSION_PATCH})
set(Boost_VERSION_STRING ${Boost_VERSION_MAJOR}_${Boost_VERSION_MINOR}_${Boost_VERSION_PATCH})
set(Boost_URL_STRING     "http://downloads.sourceforge.net/project/boost/boost/${Boost_VERSION}")
set(Boost_ARCHIVE_FILE   boost_${Boost_VERSION_STRING}.tar.gz)
set(Boost_MD5_SUM        313a11e97eb56eb7efd18325354631be) 

#
# TPL: BoostCmake
#
set(BoostCmake_VERSION_MAJOR 1)
set(BoostCmake_VERSION_MINOR 46)
set(BoostCmake_VERSION_PATCH 1)
set(BoostCmake_VERSION        ${BoostCmake_VERSION_MAJOR}.${BoostCmake_VERSION_MINOR}.${BoostCmake_VERSION_PATCH})
set(BoostCmake_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(BoostCmake_ARCHIVE_FILE   boost-cmake-cmake-${BoostCmake_VERSION}.tar.gz)
set(BoostCmake_MD5_SUM        c85221a1bf8b45537538249b9e912f4c) 

#
# TPL: HDF5
#
set(HDF5_VERSION_MAJOR 1)
set(HDF5_VERSION_MINOR 8)
set(HDF5_VERSION_PATCH 8)
set(HDF5_VERSION ${HDF5_VERSION_MAJOR}.${HDF5_VERSION_MINOR}.${HDF5_VERSION_PATCH})
set(HDF5_URL_STRING    "http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-${HDF5_VERSION}/src")
set(HDF5_ARCHIVE_FILE   hdf5-${HDF5_VERSION}.tar.gz)
set(HDF5_MD5_SUM        1196e668f5592bfb50d1de162eb16cff)      

#
# TPL: netCDF
#
set(NetCDF_VERSION_MAJOR 4)
set(NetCDF_VERSION_MINOR 1)
set(NetCDF_VERSION_PATCH 3)
set(NetCDF_VERSION ${NetCDF_VERSION_MAJOR}.${NetCDF_VERSION_MINOR}.${NetCDF_VERSION_PATCH})
set(NetCDF_URL_STRING     "http://www.unidata.ucar.edu/downloads/netcdf/ftp")
set(NetCDF_ARCHIVE_FILE   netcdf-${NetCDF_VERSION}.tar.gz)
set(NetCDF_MD5_SUM        ead16cb3b671f767396387dcb3c1a814) 

#
# TPL: ASCEMIO
#
set(ASCEMIO_VERSION_MAJOR 2)
set(ASCEMIO_VERSION_MINOR 0)
set(ASCEMIO_VERSION ${ASCEMIO_VERSION_MAJOR}.${ASCEMIO_VERSION_MINOR})
set(ASCEMIO_URL_STRING    "http://software.lanl.gov/ascem/tpls")
set(ASCEMIO_ARCHIVE_FILE   ascem-io-${ASCEMIO_VERSION}.tar.gz)
set(ASCEMIO_MD5_SUM       04d1fba6b566b38628f503a3d39c6883)      

#
# TPL: ExodusII
#
set(ExodusII_VERSION_MAJOR 5)
set(ExodusII_VERSION_MINOR 14)
set(ExodusII_VERSION ${ExodusII_VERSION_MAJOR}.${ExodusII_VERSION_MINOR})
#set(ExodusII_URL_STRING     "http://downloads.sourceforge.net/project/exodusii/exodusii/${ExodusII_VERSION}")
set(ExodusII_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(ExodusII_ARCHIVE_FILE   exodusii-${ExodusII_VERSION}.tar)
set(ExodusII_MD5_SUM        ed150bd50cbb2fb9c6d5fe95a9d46eff) 

#
# TPL: MSTK
#
set(MSTK_VERSION_MAJOR 1)
set(MSTK_VERSION_MINOR 85)
set(MSTK_VERSION_PATCH rc2)
set(MSTK_VERSION ${MSTK_VERSION_MAJOR}.${MSTK_VERSION_MINOR}${MSTK_VERSION_PATCH})
set(MSTK_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(MSTK_ARCHIVE_FILE   mstk-${MSTK_VERSION}.tgz)
set(MSTK_MD5_SUM       d2acd6da234e580b8d389224f3b2b9b1) 

#
# TPL: Trilinos
#
set(Trilinos_VERSION_MAJOR 10)
set(Trilinos_VERSION_MINOR 6)
set(Trilinos_VERSION_PATCH 4)
set(Trilinos_VERSION ${Trilinos_VERSION_MAJOR}.${Trilinos_VERSION_MINOR}.${Trilinos_VERSION_PATCH})
set(Trilinos_URL_STRING     "http://software.lanl.gov/ascem/tpls")
set(Trilinos_ARCHIVE_FILE   trilinos-${Trilinos_VERSION}-Source.tar.gz)
set(Trilinos_MD5_SUM        75b393e633bde4d9565df304f84b52e4) 
