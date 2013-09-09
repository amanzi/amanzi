if (BUILD_TESTS) 
project(ADVECTION_DIFFUSION_TEST)
    # Copy test directory files if an out of source build
    if (NOT (${ADVECTION_DIFFUSION_TEST_SOURCE_DIR} EQUAL ${ADVECTION_DIFFUSION_TEST_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy ${ADVECTION_DIFFUSION_TEST_SOURCE_DIR}/amanzi-s-1d-amr-advection.xml
               ${ADVECTION_DIFFUSION_TEST_BINARY_DIR}/amanzi-s-1d-amr-advection.xml) 
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy ${ADVECTION_DIFFUSION_TEST_SOURCE_DIR}/amanzi-s-amr-AD.xml
               ${ADVECTION_DIFFUSION_TEST_BINARY_DIR}/amanzi-s-amr-AD.xml) 
    endif()

    add_amanzi_test(examples-verification-flow-richards-1d-2a amanzi "--xml_file=amanzi-s-1d-amr-advection.xml" KIND int)
    add_amanzi_test(examples-verification-flow-richards-1d-2b amanzi "--xml_file=amanzi-s-amr-AD.xml" KIND int)

#    add_amanzi_test(examples-verification-flow-richards-1d-2a-n4 amanzi "--xml_file=case2a.xml" KIND int NPROCS 4)
#    add_amanzi_test(examples-verification-flow-richards-1d-2b-n4 amanzi "--xml_file=case2b.xml" KIND int NPROCS 4)
#    add_amanzi_test(examples-verification-flow-richards-1d-2c-n4 amanzi "--xml_file=case2c.xml" KIND int NPROCS 4)

endif()
