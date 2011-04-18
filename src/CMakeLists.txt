# -*- mode: cmake -*-

# Standard utilities -- no dependence on mesh base
add_subdirectory(error_handling)
add_subdirectory(utils)
add_subdirectory(input)

# Mesh 
add_subdirectory(mesh)
add_subdirectory(mesh_data)
add_subdirectory(simple_mesh)
add_subdirectory(exodus)

# Mesh Frameworks

# STK (Trilinos Package)
if (ENABLE_STK_Mesh)
    add_subdirectory(stk_mesh)
endif(ENABLE_STK_Mesh)    

# MOAB Mesh framework
if (ENABLE_MOAB_Mesh)
    add_subdirectory(moab_mesh)
endif(ENABLE_MOAB_Mesh)    

# MSTK Mesh framework
if (ENABLE_MSTK_Mesh)
    add_subdirectory(mstk_mesh)
endif(ENABLE_MSTK_Mesh)    

add_subdirectory(mesh_factory)
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
