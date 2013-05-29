#
# Functions for building and managing binaries.
#

include(CMakeParseArguments)
include(PrintVariable)
include(InstallManager)

#
# Usage:
#
# REGISTER_TO_FACTORY(SOURCE file1 file2 ....)
# 
#
#
# 
MACRO(REGISTER_TO_FACTORY)

  # --- Parse the input
  set(singleValueArgs SOURCE PATH)
  set(options "")
  cmake_parse_arguments(REG "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  LIST(APPEND REG_FILE_LIST ${REG_PATH}/${REG_SOURCE})

endmacro(REGISTER_TO_FACTORY)


MACRO(CREATE_REG_SOURCE)

  # --- Parse the input
  set(singleValueArgs FILENAME)
  set(options "")
  cmake_parse_arguments(REG "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})  

  foreach(_file ${REG_FILE_LIST})

    #add_custom_command(OUTPUT ${REG_FILENAME}
    #                   COMMAND /bin/cat ${_file}  >> ${REG_FILENAME}
    #                   DEPENDS ${_file}) 

    message("${_file}")

  endforeach()



endmacro(CREATE_REG_SOURCE)