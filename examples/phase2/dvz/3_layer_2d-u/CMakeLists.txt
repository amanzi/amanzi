if (BUILD_TESTS) 
project(3_LAYER_2D-U_TEST)
    # Copy test directory files if an out of source build
    if (NOT (${3_LAYER_2D-U_TEST_SOURCE_DIR} EQUAL ${3_LAYER_2D-U_TEST_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy ${3_LAYER_2D-U_TEST_SOURCE_DIR}/dvz_3_layer_2d.xml ${3_LAYER_2D-U_TEST_BINARY_DIR}/dvz_3_layer_2d.xml) 
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy ${3_LAYER_2D-U_TEST_SOURCE_DIR}/dvz_3_layer_true_2d.xml ${3_LAYER_2D-U_TEST_BINARY_DIR}/dvz_3_layer_true_2d.xml) 
    endif()

    add_amanzi_test(examples-phase2-dvz-3_layer_2d-u amanzi "--xml_file=dvz_3_layer_2d.xml" KIND int)
    SET_TESTS_PROPERTIES(examples-phase2-dvz-3_layer_2d-u PROPERTIES FAIL_REGULAR_EXPRESSION "Amanzi::SIMULATION_FAILED")

    add_amanzi_test(examples-phase2-dvz-3_layer_2d_true_2d amanzi "--xml_file=dvz_3_layer_true_2d.xml" KIND int)
    SET_TESTS_PROPERTIES(examples-phase2-dvz-3_layer_2d_true_2d PROPERTIES FAIL_REGULAR_EXPRESSION "Amanzi::SIMULATION_FAILED")
endif()
