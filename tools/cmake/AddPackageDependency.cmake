# -*- mode: cmake -*-

#
# Amanzi 
#       ADD_PACKAGE_DEPENDENCY(<PACKNAME> DEPENDS_ON <req_pack>)
#

# CMake module
include(CMakeParseArguments)

# Amanzi modules
include(PrintVariable)

function(ADD_PACKAGE_DEPENDENCY)

    # Macro: _print_usage
    macro(_print_usage)
        message("\nADD_PACKAGE_DEPENDENCY(<target_package> DEPENDS_ON <req_package>)\n"
                " Add req_package to target_package dependencies.\n")
    endmacro(_print_usage)

    # Parse the arguments
    set(_options "")
    set(_oneValue "DEPENDS_ON")
    set(_multiValue "")
    cmake_parse_arguments(ADD_PACK "${_options}" "${_oneValue}" "${_multiValue}" ${ARGN})
     
    # Define the target package name
    list(GET ADD_PACK_UNPARSED_ARGUMENTS 0 target_package)
    if ( NOT target_package )
        _print_usage()
        message(FATAL_ERROR "Must define a target_package")
    endif()

    # Define the required package
    set(req_package "")
    if(ADD_PACK_DEPENDS_ON)
        set(req_package ${ADD_PACK_DEPENDS_ON})
    else()
        _print_usage()
        message(FATAL_ERROR "Must define a required package")
    endif()    

    # Find required package
    message(STATUS "${target_package} depends on ${req_package}")
    message(STATUS "Updatig ${target_package}_LIBRARIES and ${target_package}_INCLUDE_DIRS")
    find_package(${req_package} REQUIRED)
    if( ${req_package}_LIBRARIES AND ${req_package}_INCLUDE_DIRS )
        set(_save_lib_list ${${target_package}_LIBRARIES})
        list(APPEND _save_lib_list ${${req_package}_LIBRARIES})
        set(${target_package}_LIBRARIES ${_save_lib_list} PARENT_SCOPE)
        set(_save_inc_list ${${target_package}_INCLUDE_DIRS})
        list(APPEND _save_inc_list ${${req_package}_INCLUDE_DIRS})
        set(${target_package}_INCLUDE_DIRS ${_save_inc_list} PARENT_SCOPE)
    endif()    

endfunction(ADD_PACKAGE_DEPENDENCY)
