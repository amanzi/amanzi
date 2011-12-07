# -*- mode: cmake -*-

# state
add_subdirectory(state)

# Process Kernels
add_subdirectory(pks)

# MPC
add_subdirectory(mpc)

# main (from $AMANZI_HOME/src/common/standalone_simulation_coordinator
add_subdirectory(executables)