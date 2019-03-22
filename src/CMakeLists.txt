# -*- mode: cmake -*-

# Global definitions of standards
# add_definitions(-DBOOST_NO_CXX11_VARIADIC_TEMPLATES)

# Standard utilities
add_subdirectory(error_handling)
add_subdirectory(utils)

# Functions
add_subdirectory(functions)

# Geometry
add_subdirectory(geometry)


if (ENABLE_Unstructured)
  # Mesh
  add_subdirectory(mesh)

  # Discretization
  add_subdirectory(whetstone)

  # Data structures
  add_subdirectory(data_structures)

  # Mesh functions
  add_subdirectory(mesh_functions)

  # Output
  add_subdirectory(output)

  # Solvers
  # add_subdirectory(solvers)

  # Time integration
  # add_subdirectory(time_integration)

  # Operators
  # add_subdirectory(operators)

  # Build State
  add_subdirectory(state)
endif()

