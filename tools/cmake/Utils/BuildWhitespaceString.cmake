include(PrintVariable)
include(ListJoin)

function(BUILD_WHITESPACE_STRING OUTPUT)
  list_join(VALUES ${ARGN} GLUE " " OUTPUT _TMP_STR)
  set(${OUTPUT} "${_TMP_STR}" PARENT_SCOPE)
endfunction(BUILD_WHITESPACE_STRING)  

