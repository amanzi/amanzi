if (BUILD_TESTS) 
project(8_LAYER-U_TEST)
    # Copy test directory files if an out of source build
    if (NOT (${8_LAYER-U_TEST_SOURCE_DIR} EQUAL ${8_LAYER-U_TEST_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy ${8_LAYER-U_TEST_SOURCE_DIR}/bc-cribs-1D.xml ${8_LAYER-U_TEST_BINARY_DIR}/bc-cribs-1D.xml) 
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy ${8_LAYER-U_TEST_SOURCE_DIR}/chk00030.h5 ${8_LAYER-U_TEST_BINARY_DIR}/chk00030.h5) 
    endif()

    add_amanzi_test(examples-phase2-dvz-8_layer-u amanzi "--xml_file=bc-cribs-1D.xml" KIND int)
    SET_TESTS_PROPERTIES(examples-phase2-dvz-8_layer-u PROPERTIES FAIL_REGULAR_EXPRESSION "Amanzi::SIMULATION_FAILED")

endif()
