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

# Grid-option agnostic chemistry and interfaces
add_subdirectory(common/chemistry)

if (ENABLE_ALQUIMIA)
  # Use Alquimia as a chemistry engine.
  add_definitions(-DALQUIMIA_ENABLED)
  add_subdirectory(common/alquimia)
endif()

if (ENABLE_Unstructured)
  # Mesh
  add_subdirectory(mesh)

  # Debug helpers
  add_subdirectory(debug)

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
  add_subdirectory(operators)

  # Build State
  add_subdirectory(state)
endif()

if (ENABLE_AmanziPhysicsModule)
  add_subdirectory(common/interface_platform)
endif()

if (ENABLE_Unstructured)  
  add_subdirectory(pks)

  if (ENABLE_AmanziPhysicsModule)
     # Multiprocess coordinator
     add_subdirectory(mpc)

     # unstructured grid simulation driver, and Amanzi binaries
     add_subdirectory(executables)
   endif()

   add_subdirectory(physics)
endif()

# structured grid
if (ENABLE_Structured AND ENABLE_AmanziPhysicsModule)
  add_subdirectory(structured_grid)
endif()

# standalone driver (produces mesh-independent executable)
if (ENABLE_AmanziPhysicsModule)
  add_subdirectory(common/standalone_simulation_coordinator)
endif()

# integration tests
if (ENABLE_AmanziPhysicsModule AND ENABLE_ATSPhysicsModule)
  add_subdirectory(interoperability)
endif()
