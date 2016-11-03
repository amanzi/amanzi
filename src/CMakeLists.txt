# -*- mode: cmake -*-

# operators -- layer between discretization and PK
add_subdirectory(operators)

add_subdirectory(constitutive_relations)

# Process Kernels
add_subdirectory(pks)

# top level coordinator
add_subdirectory(coordinator)

# fortran->c++ interface support
if (Amanzi_DISABLE_MPI_CXX_BINDINGS)
  add_subdirectory(xlib)
endif()


# main (from $AMANZI_HOME/src/common/standalone_simulation_coordinator
add_subdirectory(executables)
