# From AmanziConfigReport.cmake:

include(PrintVariable)

set(build_timestamp "Not available on this platform")
if (UNIX)
    execute_process(COMMAND "date"
                    RESULT_VARIABLE _ret_code
                    OUTPUT_VARIABLE _stdout
                    ERROR_VARIABLE  _stderr
                    )
    string(REGEX REPLACE "[\n\r]" "" build_timestamp ${_stdout})
endif()    


macro(_write_dir_names directories)
  foreach(dir ${${directories}})
    list(APPEND dir_flags "-L${dir} ")
  endforeach()
  file(APPEND ${AMANZI_LINK_LINE} ${dir_flags} " ")
endmacro()

function(_write_lib_names libraries)
  foreach(lib ${${libraries}})
    list(APPEND lib_names "-l${lib} ")
  endforeach()
  file(APPEND ${AMANZI_LINK_LINE} ${lib_names} " ")
endfunction()

macro(link_list_add)
  SET(package ${PROJECT_NAME})

  _write_dir_names(${package}_LIBRARY_DIR)
  _write_dir_names(${package}_LIBRARY_DIRS)

  print_variable(${package}_LIBRARIES)
  _write_lib_names(${package}_LIBRARIES)
endmacro()

macro(create_link_line)
  message(STATUS, "Writing link line to file ${AMANZI_LINK_LINE}")
  file(WRITE ${AMANZI_LINK_LINE} "-lamanzi ")
endmacro()