if (BUILD_TESTS) 
project(3_LAYER_1D-U-SSFLOW_TEST)
    # Copy test directory files if an out of source build
    if (NOT (${3_LAYER_1D-U-SSFLOW_TEST_SOURCE_DIR} EQUAL ${3_LAYER_1D-U-SSFLOW_TEST_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy ${3_LAYER_1D-U-SSFLOW_TEST_SOURCE_DIR}/dvz_3_layer_1d-ssflow.xml ${3_LAYER_1D-U-SSFLOW_TEST_BINARY_DIR}/dvz_3_layer_1d-ssflow.xml) 
    endif()

    add_amanzi_test(examples-phase2-dvz-3_layer_1d-u-ssflow amanzi "--xml_file=dvz_3_layer_1d-ssflow.xml" KIND int)
    SET_TESTS_PROPERTIES(examples-phase2-dvz-3_layer_1d-u-ssflow PROPERTIES FAIL_REGULAR_EXPRESSION "Amanzi::SIMULATION_FAILED")
endif()
