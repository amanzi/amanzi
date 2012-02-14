# -*- mode: cmake -*-

include_directories(${ATS_SOURCE_DIR}/src/data_structures)
add_library(ats_solvers NonlinarKrylovAccelerator.cc)


# if (BUILD_TESTS)
#   # Add UnitTest include directoy

#   include_directories(${Amanzi_INCLUDE_DIR})
#   include_directories(${Amanzi_TPL_Trilinos_INCLUDE_DIRS})
#   include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})
#   include_directories(${ATS_SOURCE_DIR}/src/ats_solvers)
#   include_directories(${ATS_SOURCE_DIR}/src/data_structures)
#   add_executable(test_ats_solvers test/Main.cc test/test_nka.cc)
#   target_link_libraries(test_ats_solvers ats_solvers data_structures amanzi_error_handling ${Amanzi_TPL_Trilinos_LIBRARIES} ${Amanzi_TPL_UnitTest_LIBRARIES} )
# endif()

