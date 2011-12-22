# -*- mode: cmake -*-

#
#  ATS
#    State class
#
message(STATUS "in state_dir: ATS_SOURCE_DIR=${ATS_SOURCE_DIR}")

# This state library will move to a new location
# once Markus has created the data manager.
# I use global properties since all the PROJECT_NAME
# variables in this directory are tied to the MPC.
# These properties will be replaced with the appropriate
# directory project name once it is moved.
#
# Library: state
#
#add_library(state Field.cc State.cc Vis.cc)
add_library(state Field.cc State.cc)


if (BUILD_TESTS)
  # Add UnitTest include directoy

  include_directories(${Amanzi_INCLUDE_DIR})
  include_directories(${Amanzi_TPL_Trilinos_INCLUDE_DIRS})
  include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})
  include_directories(${ATS_SOURCE_DIR}/src/state)
  add_executable(test_state test/test_copy.cc)
  target_link_libraries(test_state state amanzi_mesh amanzi_cell_geometry amanzi_mesh_factory amanzi_interface amanzi_stk_mesh amanzi_state ${Amanzi_TPL_UnitTest_LIBRARIES})
endif()