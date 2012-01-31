# -*- mode: cmake -*-

include_directories(${ATS_SOURCE_DIR}/src/data_structures)
include_directories(${ATS_SOURCE_DIR}/src/ats_solvers)
add_library(time_integration SolutionHistory.cc ImplicitTIBDF2.cc ExplicitTIRK.cc)


if (BUILD_TESTS)
  # Add UnitTest include directoy

  include_directories(${Amanzi_INCLUDE_DIR})
  include_directories(${Amanzi_TPL_Trilinos_INCLUDE_DIRS})
  include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})
  include_directories(${ATS_SOURCE_DIR}/src/time_integration)
  include_directories(${ATS_SOURCE_DIR}/src/data_structures)
  add_executable(test_time_integration test/Main.cc test/test_ode.cc test/test_explicit.cc test/test_solhist.cc)
  target_link_libraries(test_time_integration time_integration ats_solvers data_structures amanzi_error_handling ${Amanzi_TPL_Trilinos_LIBRARIES} ${Amanzi_TPL_UnitTest_LIBRARIES} )
endif()

