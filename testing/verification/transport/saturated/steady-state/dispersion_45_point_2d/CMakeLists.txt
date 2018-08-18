# -*- mode: cmake -*-

project (DISPERSION_45_POINT_2D)

set(test_name dispersion_45_point_2d)
set(compfile obs5_2_1_r5.out)
set(dataset "")
set(step 1)

if (ENABLE_TESTS)
    add_amanzi_test(run_${test_name} KIND AMANZI AMANZI_INPUT ${DISPERSION_45_POINT_2D_SOURCE_DIR}/amanzi_${test_name}-u.xml)
    add_test(NAME compare_${test_name} 
	COMMAND ${PYTHON_EXECUTABLE} ${DISPERSION_45_POINT_2D_SOURCE_DIR}/convergence_test.py ${DISPERSION_45_POINT_2D_SOURCE_DIR})
	 #COMMAND diff ${DISPERSION_45_POINT_2D_BINARY_DIR}/${compfile}
	 #                 ${DISPERSION_45_POINT_2D_SOURCE_DIR}/golden_output/${compfile} 
	 #			  #${dataset}/${step}  ${dataset}/${step} 
	 #			  )
    set_tests_properties( compare_${test_name} PROPERTIES DEPENDS run_${test_name} )
    set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "Comparison Failed" )
    #set_tests_properties( compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "differs from" )
    _append_test_label(compare_${test_name} REG)
    _append_test_label(run_${test_name} REG)
endif()
