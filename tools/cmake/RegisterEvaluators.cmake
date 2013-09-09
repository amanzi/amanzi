
# this bit of code exists to autogenerate an include file
# that contains the registration of evaluator factories
# this include file is only to be included in the file
# that contains main, and not in any file that is to be
# included in a library

function(cat IN_FILE OUT_FILE)
  file(READ ${IN_FILE} CONTENTS)
  file(APPEND ${OUT_FILE} "${CONTENTS}")
endfunction()


function(register_evaluators EVALUATOR_REG_FILE_NAME EVALUATOR_REG_LIST_NAME)

  # Prepare a temporary file to "cat" to:
  file(WRITE ${EVALUATOR_REG_FILE_NAME}.in "")

  get_property(VAR_FACTORY_REG_HEADERS GLOBAL PROPERTY EVALUATOR_REG_LIST_NAME)

  # Call the "cat" function for each input file
  foreach(HEADER ${VAR_FACTORY_REG_HEADERS})
    cat(${HEADER} ${EVALUATOR_REG_FILE_NAME}.in)
  endforeach()

  # Copy the temporary file to the final location
  configure_file(${EVALUATOR_REG_FILE_NAME}.in ${CMAKE_CURRENT_SOURCE_DIR}/${EVALUATOR_REG_FILE_NAME} COPYONLY)
  
  install(FILES ${EVALUATOR_REG_FILE_NAME}  DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

endfunction()

function(register_evaluator_with_factory EVALUATOR_STUB_NAME EVALUATOR_REG_LIST_NAME)
  get_property(VAR_TMP GLOBAL PROPERTY EVALUATOR_REG_LIST_NAME)
  set(VAR_TMP ${VAR_TMP} ${CMAKE_CURRENT_SOURCE_DIR}/${EVALUATOR_STUB_NAME})
  set_property(GLOBAL PROPERTY EVALUATOR_REG_LIST_NAME ${VAR_TMP})
endfunction()