# -*- mode: cmake -*-

# Standard utilities
add_subdirectory(error_handling)

if (ENABLE_Unstructured)
  add_subdirectory(utils)
  add_subdirectory(input)

  # Geometry
  add_subdirectory(geometry)

  # Mesh
  add_subdirectory(mesh)

  # Output
  add_subdirectory(output)

  # Solvers
  add_subdirectory(solvers)

  # Time integration
  add_subdirectory(time_integration)

  # Build State
  add_subdirectory(state)

endif()

# Grid-option agnostic chemistry and interfaces
add_subdirectory(chemistry)
add_subdirectory(common/interface_platform)

if (ENABLE_Unstructured)
  add_subdirectory(chemistry_pk)

  add_subdirectory(flow)
  add_subdirectory(transport)

  # Multiprocess coordinator
  add_subdirectory(mpc)
endif()

if (ENABLE_Unstructured)
  # unstructured grid simualtion driver, and Amanzi binaries
  add_subdirectory(executables)
endif()

# structured grid
if (ENABLE_Structured)
  add_subdirectory(structured_grid)
endif()

# standalone driver (produces mesh-independent executable)
add_subdirectory(common/standalone_simulation_coordinator)
