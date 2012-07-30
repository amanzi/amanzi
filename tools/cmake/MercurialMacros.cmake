#
# Macros for Mercurial command
#

include(CMakeParseArguments)
include(PrintVariable)

function(hg_id_command action repo output)

  find_package(Mercurial)
  if ( MERCURIAL_FOUND )
    set(cmd ${MERCURIAL_EXECUTABLE} id --${action} ${repo})
    execute_process(COMMAND ${cmd}
                    RESULT_VARIABLE err_occurred 
                    OUTPUT_VARIABLE cmd_output
                    ERROR_VARIABLE err
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_STRIP_TRAILING_WHITESPACE)
    if(err_occurred)
      message(WARNING "Could not determine mercurial ${action} in ${repo}:\n${err}")
      set(cmd_output cmd_output-NOTFOUND)
    endif()
  else()
    message(WARNING "Could not locate mercurial. Can not define mercurial ${action}")
    set(cmd_output cmd_output-NOTFOUND)
  endif()

  set(${output} ${cmd_output} PARENT_SCOPE)

endfunction(hg_id_command)


# Return branch
macro(MERCURIAL_BRANCH branch)


  set(oneValueArgs REPOSITORY)
  cmake_parse_arguments(ARGS "" "${oneValueArgs}" "" ${ARGN})

  if ( NOT ARGS_REPOSITORY )
    set(ARGS_REPOSITORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()

  hg_id_command(branch ${ARGS_REPOSITORY} ${branch})

endmacro(MERCURIAL_BRANCH)

macro(MERCURIAL_GLOBAL_ID global_id)
 
  set(oneValueArgs REPOSITORY)
  cmake_parse_arguments(ARGS "" "${oneValueArgs}" "" ${ARGN})

  if ( NOT ARGS_REPOSITORY )
    set(ARGS_REPOSITORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()

  hg_id_command(id ${ARGS_REPOSITORY} ${global_id})

endmacro(MERCURIAL_GLOBAL_ID)

macro(MERCURIAL_LOCAL_ID local_id)

  set(oneValueArgs REPOSITORY)
  cmake_parse_arguments(ARGS "" "${oneValueArgs}" "" ${ARGN})

  if ( NOT ARGS_REPOSITORY )
    set(ARGS_REPOSITORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()

  hg_id_command(num ${ARGS_REPOSITORY} ${local_id})

endmacro(MERCURIAL_LOCAL_ID)

