project(DUAL_POROSITY_1D)

set(test_name dual_porosity_1d)
set(compfile observation.out)
set(dataset "")
set(step 1)

if (ENABLE_TESTS)
  add_amanzi_test(run_${test_name} KIND AMANZI AMANZI_INPUT ${DUAL_POROSITY_1D_SOURCE_DIR}/amanzi_${test_name}-u.xml)
  add_test(NAME compare_${test_name} 
  COMMAND ${PYTHON_EXECUTABLE} ${DUAL_POROSITY_1D_SOURCE_DIR}/convergence_test.py ${DUAL_POROSITY_1D_SOURCE_DIR})
  #COMMAND diff ${DUAL_POROSITY_1D_BINARY_DIR}/${compfile}
  #             ${DUAL_POROSITY_1D_SOURCE_DIR}/golden_output/${compfile} 
  #	        #${dataset}/${step}  ${dataset}/${step} 
  set_tests_properties(compare_${test_name} PROPERTIES DEPENDS run_${test_name})
  set_tests_properties(compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "Comparison Failed")
  #set_tests_properties(compare_${test_name} PROPERTIES FAIL_REGULAR_EXPRSSION "differs from")
  _append_test_label(compare_${test_name} REG)
  _append_test_label(run_${test_name} REG)
endif()
