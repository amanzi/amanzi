#
# CHECK_DOWNLOAD
#
# USAGE:
#  CHECK_DOWNLOAD(TEST_URL web-address
#                 TEST_FILE file
#                 TIMEOUT seconds)
#
# Attempt to download file from web-address using the system curl binary.
# Default TIMEOUT for this command is 60 seconds. If download fails, throws
# fatal error. 
#
include(CMakeParseArguments)
include(PrintVariable)

function(CHECK_DOWNLOAD)


  # Parse the arguments
  set(_flags    "")
  set(_oneValue "TEST_URL;TEST_FILE;TIMEOUT")
  set(_multiValue "")
  cmake_parse_arguments(PARSE "${_flags}" "${_oneValue}" "${_multiValue}" ${ARGN})

  # Default timeout is 60 seconds
  set(command_timeout 60)
  if ( PARSE_TIMEOUT )
    set(command_timeout ${PARSE_TIMEOUT})
  endif()  

  # Need a URL AND FILE name
  if ( PARSE_TEST_URL AND PARSE_TEST_FILE )
    set(url_string ${PARSE_TEST_URL}/${PARSE_TEST_FILE})
  else()
    message(FATAL_ERROR "Invalid arguments to CHECK_DOWNLOAD. "
                        "Must define TEST_URL AND TEST_FILE")
  endif()                    
       
  message(STATUS "Checking external downloads (${url_string})")

  file(DOWNLOAD
       ${url_string}
       ${CMAKE_CURRENT_BINARY_DIR}/${PARSE_TEST_FILE}
    TIMEOUT ${command_timeout}
       STATUS result)
  list(GET result 0 ret_code)
  list(GET result 1 error_str)
  if ( "${ret_code}" EQUAL 0 )
    message(STATUS "Checking external downloads (${url_string}) -- works")
  else()  
    message(SEND_ERROR "Failed to download ${url_string} "
                       "Return Code:${ret_code}\n"
                       "Output:\n${error_str}\n"
                       "You can disable external downloads with"
                       "\n-DDISABLE_EXTERNAL_DOWNLOAD:BOOL=TRUE\n"
                       "If external downloads are disabled, ALL TPL source files must be found in one directory"
                       " and define this directory with"
                       "\n-D TPL_DOWNLOAD_DIR:FILEPATH=\n")
    message(FATAL_ERROR "Failed to download ${PARSE_TEST_FILE} from ${PARSE_TEST_URL}")

  endif()  


endfunction(CHECK_DOWNLOAD)
