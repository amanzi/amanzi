# -*- mode: cmake -*-

# Standard utilities -- no dependence on mesh base
add_subdirectory(error_handling)
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

add_subdirectory(chemistry)
add_subdirectory(flow)
add_subdirectory(transport)

# Multiprocess coordinator
add_subdirectory(mpc)

# Amanzi binaries
add_subdirectory(executables)
