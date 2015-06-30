# -*- mode: cmake -*-

project (BUTLER_POD_2D)

set(test_name butler_pod_2d)
set(compfile observations.out)
set(dataset "")
set(step 1)

if (ENABLE_TESTS)
    add_amanzi_comparison_test(${test_name} 
                               INPUT amanzi_${test_name}.xml
                               FILES mesh_cylinder.exo
                               TOLERANCE 2e-3
                               OBSERVATION Drawdown
                               REFERENCE ${BUTLER_POD_2D_SOURCE_DIR}/golden_output/${compfile})
endif()
#    add_amanzi_test(run_${test_name} KIND AMANZI AMANZI_INPUT ${BUTLER_POD_2D_SOURCE_DIR}/amanzi_${test_name}.xml)
#    add_test(NAME compare_${test_name} 
#             COMMAND ${PYTHON_EXECUTABLE} ${BUTLER_POD_2D_SOURCE_DIR}/convergence_test.py ${BUTLER_POD_2D_SOURCE_DIR})
#	 #COMMAND diff ${BUTLER_POD_2D_BINARY_DIR}/${compfile}
#	 #                 ${BUTLER_POD_2D_SOURCE_DIR}/golden_output/${compfile} 
#	 #			  #${dataset}/${step}  ${dataset}/${step} 
#	 #			  )
#    set_tests_properties( compare_${test_name} PROPERTIES DEPENDS run_${test_name} )
#    set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "Comparison Failed" )
#    #set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "differs from" )
#    _append_test_label(compare_${test_name} REG)
#    _append_test_label(run_${test_name} REG)
#endif()
