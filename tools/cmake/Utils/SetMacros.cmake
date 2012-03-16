#
# SetMacros
#
#  Collection of useful macros that append variables and 
#  define global variables.


#
# GLOBAL_SET(VAR_NAME ... )
#  Set global variable VAR_NAME
#
macro (GLOBAL_SET VAR)
  set(${VAR} ${ARGN} CACHE INTERNAL "")
endmacro(GLOBAL_SET) 

#
# APPEND_SET(VAR item1 [item2] ...)
#  Append item1 [item2] ... to variable VAR
# 
macro(APPEND_SET VAR)
  set(${VAR} ${${VAR}} ${ARGN})
endmacro(APPEND_SET)

#
# PREPEND_SET(VAR item1 [item2] ...)
#   Prepend item1 [item2] ... to variable VAR
#
macro(PREPEND_SET VAR)
  set(${VAR} ${ARGN} ${${VAR}})
endmacro(PREPEND_SET)  

#
# GLOBAL_APPEND_SET(VAR item1 [item2] ...)
#  Append item1 [item2] ... to global variable VAR
# 
macro(GLOBAL_APPEND_SET VAR)
  global_set(${VAR} ${${VAR}} ${ARGN})
endmacro(GLOBAL_APPEND_SET)

#
# GLOBAL_PREPEND_SET(VAR item1 [item2] ...) 
#  Prepend item1 [item2] ... to global variable VAR
#
macro(GLOBAL_PREPEND_SET VAR)
  global_set(${VAR} ${ARGN} ${${VAR}})
endmacro(GLOBAL_PREPEND_SET)  

