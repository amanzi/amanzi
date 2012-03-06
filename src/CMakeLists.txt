# -*- mode: cmake -*-

# state
add_subdirectory(data_structures)
add_subdirectory(state)

# Process Kernels
add_subdirectory(pks)

# top level coordinator
add_subdirectory(coordinator)

# solvers with ATS data structures
add_subdirectory(ats_solvers)

# time integration with ATS data structures
add_subdirectory(time_integration)

# main (from $AMANZI_HOME/src/common/standalone_simulation_coordinator
add_subdirectory(executables)
