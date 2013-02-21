if (BUILD_TESTS) 
project(RICHARDS_1D_TEST)
    # Copy test directory files if an out of source build
    if (NOT (${RICHARDS_1D_TEST_SOURCE_DIR} EQUAL ${RICHARDS_1D_TEST_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy ${RICHARDS_1D_TEST_SOURCE_DIR}/case2a.xml ${RICHARDS_1D_TEST_BINARY_DIR}/case2a.xml) 
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy ${RICHARDS_1D_TEST_SOURCE_DIR}/case2b.xml ${RICHARDS_1D_TEST_BINARY_DIR}/case2b.xml) 
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy ${RICHARDS_1D_TEST_SOURCE_DIR}/case2c.xml ${RICHARDS_1D_TEST_BINARY_DIR}/case2c.xml) 
    endif()

    add_amanzi_test(examples-verification-flow-richards-1d-2a amanzi "--xml_file=case2a.xml" KIND int)
    add_amanzi_test(examples-verification-flow-richards-1d-2b amanzi "--xml_file=case2b.xml" KIND int)
    add_amanzi_test(examples-verification-flow-richards-1d-2c amanzi "--xml_file=case2c.xml" KIND int)

#    add_amanzi_test(examples-verification-flow-richards-1d-2a-n2 amanzi "--xml_file=case2a.xml" KIND int NPROCS 2)
#    add_amanzi_test(examples-verification-flow-richards-1d-2b-n2 amanzi "--xml_file=case2b.xml" KIND int NPROCS 2)
#    add_amanzi_test(examples-verification-flow-richards-1d-2c-n2 amanzi "--xml_file=case2c.xml" KIND int NPROCS 2)

#    add_amanzi_test(examples-verification-flow-richards-1d-2a-n4 amanzi "--xml_file=case2a.xml" KIND int NPROCS 4)
#    add_amanzi_test(examples-verification-flow-richards-1d-2b-n4 amanzi "--xml_file=case2b.xml" KIND int NPROCS 4)
#    add_amanzi_test(examples-verification-flow-richards-1d-2c-n4 amanzi "--xml_file=case2c.xml" KIND int NPROCS 4)

endif()
