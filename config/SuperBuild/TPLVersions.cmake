#  -*- mode: cmake -*-

#
# TPLVersions
#    Define the versions, approved download locations for each TPL
#
#


#
# TPL: CURL
#
set(CURL_VERSION_MAJOR 7)
set(CURL_VERSION_MINOR 21)
set(CURL_VERSION_PATCH 6)
set(CURL_VERSION ${CURL_VERSION_MAJOR}.${CURL_VERSION_MINOR}.${CURL_VERSION_PATCH})
set(CURL_URL_STRING     "http://curl.haxx.se/download")
set(CURL_ARCHIVE_FILE   curl-${CURL_VERSION}.tar.gz)
set(CURL_MD5_SUM       c502b67898b4a1bd687fe1b86419a44b) 

if(TPL_DOWNLOAD_DIR)
  set(CURL_DOWNLOAD_LOCATION ${TPL_DOWNLOAD_DIR} )
else()  
  set(CURL_DOWNLOAD_LOCATION ${CURL_URL_STRING})
endif()  
