# -*- mode: cmake -*-

#
# Amanzi 
#       Python CMakeLists.txt file
#

project(AmanziPython)

# Amanzi CMake modules
include(PrintVariable)
include(InstallManager)

# Search for Python
find_package(Python)

if ( PYTHON_FOUND )

  # Define the install location
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}')"
                  OUTPUT_VARIABLE AmanziPython_INSTALL_PREFIX
                  RESULT_VARIABLE prefix_err
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(prefix_err)
    message(SEND_ERROR "Failed to define the python from distutils. Set to ${CMAKE_INSTALL_PREFIX}/python.")
    set(AmanziPython_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/python)
  endif()

  # Define the Amanzi binary install location
  if ( TARGET amanzi )
    get_target_property(base amanzi OUTPUT_NAME)

    set(Amanzi_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/bin/${base})

  endif()

  # List of scripts to install
  set(python_scripts)

  # Configure files
  configure_file(${AmanziPython_SOURCE_DIR}/vv_driver.py
                 ${AmanziPython_BINARY_DIR}/vv_driver.py
                 @ONLY)

  list(APPEND python_scripts ${AmanziPython_BINARY_DIR}/vv_driver.py)             

               


  # List of directories to install
  set(python_module_directories amanzi)

  message(STATUS "Install Python tools and modules in: ${AmanziPython_INSTALL_PREFIX}")

  # Install the script in the python library for now
  INSTALL(FILES ${python_scripts} DESTINATION ${AmanziPython_INSTALL_PREFIX}
          PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                      GROUP_EXECUTE GROUP_READ)

  # Install the module directories
  foreach ( dir ${python_module_directories} )
    install(DIRECTORY ${dir} DESTINATION ${AmanziPython_INSTALL_PREFIX}
            PATTERN "*.pyc" EXCLUDE
            PATTERN "${dir}/*.py"
            PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                        GROUP_EXECUTE GROUP_READ)
  endforeach()



endif()


