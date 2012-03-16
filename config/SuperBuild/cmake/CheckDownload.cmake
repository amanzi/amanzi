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

  # Use the system curl executable
  find_program(_CURL_EXECUTABLE curl)
  
  message(STATUS "Checking external downloads with ${_CURL_EXECUTABLE}")

  # Parse the arguments
  set(_flags    "")
  set(_oneValue "TEST_URL;TEST_FILE;TIMEOUT")
  set(_multiValue "")
  cmake_parse_arguments(PARSE "${_flags}" "${_oneValue}" "${_multiValue}" ${ARGN})

  set(command_timeout 60)
  if ( PARSE_TIMEOUT )
    set(command_timeout ${PARSE_TIMEOUT})
  endif()  

  execute_process(COMMAND ${_CURL_EXECUTABLE} ${PARSE_TEST_URL}/${PARSE_TEST_FILE}
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    TIMEOUT ${command_timeout}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE command_output
    ERROR_VARIABLE command_output
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
    )

  if ( NOT "${result}" EQUAL 0 )
    message(SEND_ERROR "Download command ${_CURL_EXECUTABLE} ${PARSE_TEST_URL}/${PARSE_TEST_FILE} failed"
                       "Output:\n${command_output} ${result}\n"
                       "You can disable external downloads with"
                       "\n-DDISABLE_EXTERNAL_DOWNLOAD:BOOL=TRUE\n"
                       "If external downloads are disabled, ALL TPL source files must be found in one directory"
                       " and define this directory with"
                       "\n-D TPL_DOWNLOAD_DIR:FILEPATH=\n")
    message(FATAL_ERROR "Failed to download ${PARSE_TEST_FILE} from ${PARSE_TEST_URL}")


  else()
    message(STATUS "Checking external downloads with ${_CURL_EXECUTABLE} -- works")
  endif()  

  set(_CURL_EXECUTABLE)        

endfunction(CHECK_DOWNLOAD)
