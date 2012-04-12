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
include_directories(${ATS_SOURCE_DIR}/src/data_structures)
include_directories(${ATS_SOURCE_DIR}/src/state)
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy)
include_directories(${ATS_SOURCE_DIR}/src/coordinator)

add_executable(ats amanzi_unstructured_grid_simulation_driver.cc main.cc )

set(pk_factory -Wl,--whole-archive pk_factory -Wl,--no-whole-archive)

target_link_libraries(ats
                      coordinator
                      ${pk_factory}
                      mpc
                      amanzi_solvers
                      amanzi_time_integration
                      amanzi_mesh
                      amanzi_cell_geometry
                      amanzi_mesh_factory
                      amanzi_input
                      amanzi_output
                      amanzi_interface
                      amanzi_stk_mesh
                      amanzi_state
                      )

install(TARGETS ats DESTINATION ${ATS_SOURCE_DIR}/bin)






