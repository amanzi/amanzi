# ############################################################################ #
#
# CMake FindPython module
#
# ############################################################################ #

include(FindPackageHandleStandardArgs)

# Search for the Python executable
if ( NOT PYTHON_EXECUTABLE )

  # Call PythonInterp to find the python executable
  find_package(PythonInterp)
  print_variable(PYTHON_EXECTUABLE)
  print_variable(PYTHON_VERSION_STRING)

endif()

# Define the version
if (PYTHON_EXECUTABLE AND (NOT PYTHON_VERSION_STRING) )
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print('%s.%s.%s' % sys.version_info[0:3])"
                  OUTPUT_VARIABLE PYTHON_VERSION_STRING
                  RESULT_VARIABLE ret
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  
  if(ret)
    message(SEND_ERROR "Failed to define PYTHON_VERSION_STRING")
  endif()

  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print('%s' % sys.version_info[0])"
    OUTPUT_VARIABLE PYTHON_VERSION_MAJOR
                  RESULT_VARIABLE ret
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  
  if(ret)
    message(SEND_ERROR "Failed to define PYTHON_VERSION_MAJOR")
  endif()

  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print('%s' % sys.version_info[1])"
    OUTPUT_VARIABLE PYTHON_VERSION_MINOR
                  RESULT_VARIABLE ret
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  
  if(ret)
    message(SEND_ERROR "Failed to define PYTHON_VERSION_MINOR")
  endif()

  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print('%s' % sys.version_info[2])"
    OUTPUT_VARIABLE PYTHON_VERSION_PATCH
                  RESULT_VARIABLE ret
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  
  if(ret)
    message(SEND_ERROR "Failed to define PYTHON_VERSION_PATCH")
  endif()


endif()

# Search for the PYTHON_INCLUDE_DIRS and PYTHON_LIBRARIES
if ( PYTHON_EXECUTABLE AND (NOT PYTHON_INCLUDE_DIRS) )

  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print sys.prefix"
                  OUTPUT_VARIABLE PYTHON_PREFIX
                  RESULT_VARIABLE ret
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  
  if ( ret )
    message(SEND_ERROR "Failed to locate Python install prefix")
  endif()

  if(PYTHON_PREFIX)
    set(_python_search_paths ${PYTHON_PREFIX}/include)
    set(_python_suffixes  ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})

    find_path(PYTHON_INCLUDE_DIRS
              Python.h
              PATHS ${PYTHON_PREFIX}/include
              PATH_SUFFIXES python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}
              NO_DEFAULT_PATH)
  endif()
                
                  



endif()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(Python DEFAULT_MSG 
                                  PYTHON_EXECUTABLE)



