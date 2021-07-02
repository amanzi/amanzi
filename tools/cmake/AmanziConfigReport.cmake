# -*- mode: cmake -*-

# #############################################################################
#    
# Amanzi Configuration Report
#
# #############################################################################
include(FeatureSummary)


# Grab global information about enabled languages and packages
# that have returned successfully from a find_package call
get_property(lang_enabled GLOBAL PROPERTY ENABLED_LANGUAGES)
get_property(pack_found GLOBAL PROPERTY PACKAGES_FOUND)
get_property(pack_not_found GLOBAL PROPERTY PACKAGES_NOT_FOUND)

# Define the build timestamp
set(build_timestamp "Not available on this platform")
if (UNIX)
    execute_process(COMMAND "date"
                    RESULT_VARIABLE _ret_code
                    OUTPUT_VARIABLE _stdout
                    ERROR_VARIABLE  _stderr
                    )
    string(REGEX REPLACE "[\n\r]" "" build_timestamp ${_stdout})
endif()    

# Useful macros
macro(_write_to_log line)
    file(APPEND ${AMANZI_CONFIG_LOG} "${line}\n")
endmacro(_write_to_log)

macro(_write_blank_line)
    _write_to_log("")
endmacro(_write_blank_line)    

macro(_create_log_file timestamp)
    file(WRITE ${AMANZI_CONFIG_LOG} "${PROJECT_NAME} Configuration\n")
    _write_to_log("Version ${AMANZI_VERSION}")
    _write_blank_line()
    _write_to_log("Build timestamp: ${timestamp}")
    _write_blank_line()
endmacro(_create_log_file)    

# Create the log file
_create_log_file(${build_timestamp})

# Write System Information
_write_blank_line()
_write_to_log( "System Information")
_write_to_log( "\tSystem:         ${CMAKE_SYSTEM}")
_write_to_log( "\tSystem Name:    ${CMAKE_SYSTEM_NAME}")
_write_to_log( "\tSystem Version: ${CMAKE_SYSTEM_VERSION}")
if ( APPLE )
    _write_to_log( "\tSystem Type:    Mac OSX")
endif(APPLE)    
if ( WIN32 )
    _write_to_log( "\tSystem Type:    Windows")
endif(WIN32)    
if ( UNIX )
    _write_to_log( "\tSystem Type:    UNIX based platform")
endif(UNIX)   

# Build directory
_write_blank_line()
_write_to_log("Build directory:\t\t${AMANZI_BINARY_DIR}")
_write_blank_line()

# Write Compiler Information
_write_blank_line()
_write_to_log( "Compilers")
_write_to_log( "\tEnabled Languages: ${lang_enabled}")
_write_to_log( "\tC COMPILER:\t\t${CMAKE_C_COMPILER}")
_write_to_log( "\tC COMPILER ID:\t\t${CMAKE_C_COMPILER_ID}")
_write_to_log( "\tCXX COMPILER:\t\t${CMAKE_CXX_COMPILER}")
_write_to_log( "\tCXX COMPILER ID:\t\t${CMAKE_CXX_COMPILER_ID}")
if (CMAKE_Fortran_COMPILER)
    _write_to_log( "\tFortran Compiler:\t\t${CMAKE_Fortran_COMPILER}")
    _write_to_log( "\tFortran Compiler ID:\t\t${CMAKE_Fortran_COMPILER_ID}")
endif()    
_write_blank_line()
set(cmake_build_type_uc)
if(CMAKE_BUILD_TYPE)
  _write_to_log( "Build type: ${CMAKE_BUILD_TYPE}")
  string(TOUPPER "${CMAKE_BUILD_TYPE}" cmake_build_type_uc)
else()
  _write_to_log("CMAKE_BUILD_TYPE was not set")
endif()
_write_to_log( "")
_write_to_log( "Compile Flags")
_write_to_log("\tC Flags")
_write_to_log("\t\tCMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
if(cmake_build_type_uc)
  _write_to_log("\t\tCMAKE_C_FLAGS_${cmake_build_type_uc}: ${CMAKE_C_FLAGS_${cmake_build_type_uc}}")
endif()
_write_to_log("\tCXX Flags")
_write_to_log("\t\tCMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
if(cmake_build_type_uc)
  _write_to_log("\t\tCMAKE_CXX_FLAGS_${cmake_build_type_uc}: ${CMAKE_CXX_FLAGS_${cmake_build_type_uc}}")
endif()
if(CMAKE_Fortran_COMPILER)
  _write_to_log("\tFortran Flags")
  _write_to_log("\t\tCMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
  if(cmake_build_type_uc)
    _write_to_log("\t\tCMAKE_Fortran_FLAGS_${cmake_build_type_uc}: ${CMAKE_Fortran_FLAGS_${cmake_build_type_uc}}")
  endif()
endif()
_write_blank_line()

_write_blank_line()
_write_to_log("Repository Information")
_write_to_log("\tSource Directory:    ${AMANZI_SOURCE_DIR}")
_write_to_log("\tGit Branch:          ${AMANZI_GIT_BRANCH}")
_write_to_log("\tGit Global ID:       ${AMANZI_GIT_GLOBAL_HASH}")
_write_blank_line()

# Write Package Information
_write_to_log("Third Party Libraries")
if(pack_found)
    list(LENGTH pack_found num_pack_found)
    _write_to_log("\tNumber Of Packages Found: ${num_pack_found}")
    _write_to_log("\tPackages: ${pack_found}")
    _write_blank_line()
    set(package_vars VERSION VERSION_STRING ROOT DIR INCLUDE_DIRS INCLUDE_DIRS LIBRARY LIBRARIES LIBRARY_DIR EXECUTABLE)
    foreach(pack ${pack_found})
        _write_to_log( "\t${pack}")
        foreach(append_var ${package_vars})
          if (${pack}_${append_var})
            _write_to_log("\t\t${pack}_${append_var}=${${pack}_${append_var}}")
          endif()
        endforeach()
        _write_blank_line()
    endforeach()    
    _write_blank_line()
endif(pack_found)
if(pack_not_found)
    _write_blank_line()
    list(LENGTH pack_not_found num_pack_not_found)
    _write_to_log("\tNumber Of Packages Not Found: ${num_pack_not_found}") 
    if ( num_pack_not_found GREATER 0 )
        _write_to_log("\tDid not find packages: ${pack_not_found}")
    else()
        _write_to_log("\tLocated all required and requested packages")
    endif()    
else()    
    _write_blank_line()
    _write_to_log("\tLocated all required and requested packages")
endif()

# Write Mesh Framework Information
_write_blank_line()
_write_to_log("Enabled mesh frameworks")
_write_to_log("\tENABLE_MESH_MOAB:   ${ENABLE_MESH_MOAB}")
_write_to_log("\tENABLE_MESH_MSTK:   ${ENABLE_MESH_MSTK}")
_write_blank_line()
_write_to_log("Enabled mesh capability")
_write_to_log("\tENABLE_Unstructured:    ${ENABLE_Unstructured}")
_write_to_log("\tENABLE_Structured:      ${ENABLE_Structured}")
_write_blank_line()

# Write CCSE configuration info
_write_to_log( "CCSE configuration")
_write_to_log( "\tENABLE_MPI:      ${ENABLE_MPI}")
_write_to_log( "\tENABLE_OpenMP:   ${ENABLE_OpenMP}")
_write_to_log( "\tSPACEDIM:        ${AMANZI_SPACEDIM}")
_write_to_log( "\tPRECISION:       ${AMANZI_PRECISION}")
_write_blank_line()

# Write install target information
_write_to_log("Installation paths")
_write_to_log("\tAMANZI_INSTALL_PREFIX       ${CMAKE_INSTALL_PREFIX}")
_write_to_log("\tAMANZI_INSTALL_INCLUDE_PATH ${CMAKE_INSTALL_PREFIX}/include")
_write_to_log("\tAMANZI_INSTALL_LIB_PATH     ${CMAKE_INSTALL_PREFIX}/lib")
_write_to_log("\tAMANZI_INSTALL_BIN_PATH     ${CMAKE_INSTALL_PREFIX}/bin")
_write_blank_line()




#print_enabled_features()
#print_disabled_features()


