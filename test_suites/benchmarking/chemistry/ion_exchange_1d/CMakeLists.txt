if (BUILD_TESTS) 
project(CHEMISTRY_1D_ION_EXCHANGE_TEST)
    # Copy test directory files if an out of source build
    if (NOT (${CHEMISTRY_1D_ION_EXCHANGE_TEST_SOURCE_DIR} EQUAL ${CHEMISTRY_1D_ION_EXCHANGE_TEST_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy ${CHEMISTRY_1D_ION_EXCHANGE_TEST_SOURCE_DIR}/amanzi-u-1d-ion-exchange.xml ${CHEMISTRY_1D_ION_EXCHANGE_TEST_BINARY_DIR}/amanzi-u-1d-ion-exchange.xml) 
    endif()

    add_amanzi_test(examples-phase2-chemistry-1d-ion-exchange amanzi "--xml_file=amanzi-u-1d-ion-exchange.xml" KIND int)
    SET_TESTS_PROPERTIES(examples-phase2-chemistry-1d-ion-exchange PROPERTIES FAIL_REGULAR_EXPRESSION "Amanzi::SIMULATION_FAILED")
endif()
