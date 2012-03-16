#  -*- mode: cmake -*-

#
# Amanzi SuperBuild
#       Amanzi TPL Build Targets
#
project(SuperBuildDoc)

# Build the docs
find_package(LATEX)


if (LATEX_COMPILER)

  if ( NOT ( ${SuperBuildDoc_SOURCE_DIR} EQUAL ${SuperBuildDoc_BINARY_DIR}) )
    execute_process(COMMAND ${CMAKE_COMMAND} -E
            copy_directory ${SuperBuildDoc_SOURCE_DIR}/figures ${SuperBuildDoc_BINARY_DIR}/figures
	    )
  endif()	  

  ADD_CUSTOM_COMMAND(
    OUTPUT        ${SuperBuildDoc_BINARY_DIR}/build_system_requirements.pdf
    COMMAND       ${PDFLATEX_COMPILER}
    ARGS          ${SuperBuildDoc_SOURCE_DIR}/build_system_requirements.tex
    DEPENDS       ${SuperBuildDoc_SOURCE_DIR}/build_system_requirements.tex
    COMMENT       "Build PDF requirements document"
    )

  ADD_CUSTOM_TARGET(LatexDoc ALL echo
    DEPENDS ${SuperBuildDoc_BINARY_DIR}/build_system_requirements.pdf
  )
else()
  message(WARNING "Could not locate LATEX compiler. Will not build"
                  " documentation")
endif()







 
