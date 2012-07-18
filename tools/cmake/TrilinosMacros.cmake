# ############################################################################ #
#
# Useful macros to be used with the Trilinos software package
#
# ############################################################################ #

#
# This macro searches for a Trilnos package if
# the ${package}_FOUND is not set. Loops through
# the ${package}_TPL_LIST and creates 
# ${package}_ENABLE_<tpl> variables
#
macro(TRILINOS_PACKAGE_ENABLED_TPLS package)

  if (NOT ${${package}_FOUND} )
    find_package(${package}
                 NO_MODULE
                 HINTS ${Trilinos_DIR}
                 PATH_SUFFIXES include)
  endif() 

  if ( ${${package}_FOUND} ) 

    foreach ( tpl ${${package}_TPL_LIST} )
      set(${package}_ENABLE_${tpl}  TRUE)
    endforeach()

  else()
    message(WARNING "Could not locate Trilinos package ${package}")
  endif()  

endmacro(TRILINOS_PACKAGE_ENABLED_TPLS)
