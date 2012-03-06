# -*- mode: cmake -*-

#
#  ATS
#    MPC class
#

# This state library will move to a new location
# once Markus has created the data manager.
# I use global properties since all the PROJECT_NAME
# variables in this directory are tied to the MPC.
# These properties will be replaced with the appropriate
# directory project name once it is moved.
#
# Library: state
#
link_directories(${ATS_BINARY_DIR}/data_structures)
link_directories(${ATS_BINARY_DIR}/state)
link_directories(${ATS_BINARY_DIR}/pks)
link_directories(${ATS_BINARY_DIR}/pks/mpc)
link_directories(${ATS_BINARY_DIR}/pks/energy)
link_directories(${ATS_BINARY_DIR}/coordinator)

include_directories(${ATS_SOURCE_DIR}/src/data_structures)
include_directories(${ATS_SOURCE_DIR}/src/state)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy)
include_directories(${ATS_SOURCE_DIR}/src/coordinator)

add_executable(ats main.cc amanzi_unstructured_grid_simulation_driver.cc)

target_link_libraries(ats
                      coordinator
                      pks
                      mpc
                      energy
                      state
                      time_integration
                      ats_solvers
                      data_structures
                      amanzi_solvers
                      amanzi_time_integration
                      amanzi_mesh
                      amanzi_cell_geometry
                      amanzi_mesh_factory
                      amanzi_flow_implementation
                      amanzi_input
                      amanzi_output
                      amanzi_interface
                      amanzi_stk_mesh
                      amanzi_state
                      )

install(TARGETS ats DESTINATION ${ATS_SOURCE_DIR}/bin)






