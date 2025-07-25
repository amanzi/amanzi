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
  get_filename_component(real_install_prefix ${CMAKE_INSTALL_PREFIX} REALPATH)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_path('purelib', vars={'base': '${real_install_prefix}'}))"
                  OUTPUT_VARIABLE AmanziPython_INSTALL_PREFIX
                  RESULT_VARIABLE prefix_err
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  
  if(prefix_err)
    message(SEND_ERROR "Failed to define the python from sysconfig. Set to ${real_install_prefix}/python.")
    set(AmanziPython_INSTALL_PREFIX ${real_install_prefix}/python)
  endif()

  # Define the Amanzi binary install location
  if ( TARGET amanzi )
    get_target_property(base amanzi OUTPUT_NAME)
    set(Amanzi_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/bin/${base})
  endif()

  if ( TARGET unscramble_viz )
    get_target_property(base unscramble_viz OUTPUT_NAME)
    set(Amanzi_UnscrambleViz_BINARY ${CMAKE_INSTALL_PREFIX}/bin/${base})
  endif()

  if ( TARGET unscramble_restart )
    get_target_property(base unscramble_restart OUTPUT_NAME)
    set(Amanzi_UnscrambleRestart_BINARY ${CMAKE_INSTALL_PREFIX}/bin/${base})
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


