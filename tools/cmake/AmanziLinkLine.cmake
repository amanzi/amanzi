include(ParseLibraryList)
include(PrintVariable)

# From AmanziConfigReport.cmake:
set(build_timestamp "Not available on this platform")
if (UNIX)
  execute_process(COMMAND "date"
    RESULT_VARIABLE _ret_code
    OUTPUT_VARIABLE _stdout
    ERROR_VARIABLE  _stderr
    )
  string(REGEX REPLACE "[\n\r]" "" build_timestamp ${_stdout})
endif()

macro(_add_item value)
  list(APPEND value_list "${value} ")
endmacro()

macro(_add_directories directories)
  foreach(directory ${directories})
    _add_item("-L${directory}")
  endforeach()
endmacro()

macro(_add_library library)
  if(EXISTS ${library})
    _add_item("${library}")
  else()
    _add_item("-l${library}")
  endif()
endmacro()


macro(_add_libraries libraries)
foreach(library ${libraries})
  _add_library(${library})
endforeach()
endmacro()



function(_parse_add_libraries library_list libraries_to_add)

if (library_list)
  parse_library_list(
    ${library_list}
    FOUND   libraries_split
    DEBUG   debug_libraries
    OPT     opt_libraries
    GENERAL general_libraries)

  if (libraries_split)
    message("Libraries for ${package} were split")
    if(${CMAKE_BUILD_TYPE} MATCHES "debug")
      message("Adding debug libraries")
      set(${libraries_to_add} "${debug_libraries}" PARENT_SCOPE)
    else()
      message("Adding optimized libraries")
      set(${libraries_to_add} "${opt_libraries}" PARENT_SCOPE)
    endif()
  else()
    set(${libraries_to_add} "${library_list}" PARENT_SCOPE)
  endif()

endif()
endfunction()



macro(link_list_add parent_library)
  SET(package ${PROJECT_NAME})

  _add_item("-l${parent_library}")

  _add_directories("${${package}_LIBRARY_DIR}")
  _add_directories("${${package}_LIBRARY_DIRS}")

  _parse_add_libraries("${${package}_LIBRARIES}" to_add)
  _add_libraries("${to_add}")

  file(APPEND ${AMANZI_LINK_LINE} ${value_list} " ")

endmacro()

macro(create_link_line)
  message(STATUS "Writing link line to file ${AMANZI_LINK_LINE}")
  file(WRITE ${AMANZI_LINK_LINE})
endmacro()