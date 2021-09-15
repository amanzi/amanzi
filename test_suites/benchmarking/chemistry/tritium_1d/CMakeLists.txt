if (BUILD_TESTS) 
project(CHEMISTRY_1D_TRITIUM_TEST)
    # Copy test directory files if an out of source build
    if (NOT (${CHEMISTRY_1D_TRITIUM_TEST_SOURCE_DIR} EQUAL ${CHEMISTRY_1D_TRITIUM_TEST_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy ${CHEMISTRY_1D_TRITIUM_TEST_SOURCE_DIR}/amanzi-u-1d-tritium.xml ${CHEMISTRY_1D_TRITIUM_TEST_BINARY_DIR}/amanzi-u-1d-tritium.xml) 
    endif()

    add_amanzi_test(examples-phase2-chemistry-1d-tritium amanzi "--xml_file=amanzi-u-1d-tritium.xml" KIND int)
    SET_TESTS_PROPERTIES(examples-phase2-chemistry-1d-tritium PROPERTIES FAIL_REGULAR_EXPRESSION "Amanzi::SIMULATION_FAILED")
endif()
