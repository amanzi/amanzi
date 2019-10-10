# -*- mode: cmake -*-

# Standard utilities
add_subdirectory(error_handling)
add_subdirectory(utils)

# Functions
add_subdirectory(functions)

# Geometry
add_subdirectory(geometry)

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
add_subdirectory(solvers)

# Time integration
add_subdirectory(time_integration)

# Operators
# add_subdirectory(operators)

# Build State
add_subdirectory(state)

# PKs
#add_subdirectory(pks)


