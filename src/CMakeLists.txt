# -*- mode: cmake -*-

# operators -- layer between discretization and PK
add_subdirectory(operators)

# Process Kernels
add_subdirectory(pks)

# top level coordinator
add_subdirectory(coordinator)

# main (from $AMANZI_HOME/src/common/standalone_simulation_coordinator
add_subdirectory(executables)
