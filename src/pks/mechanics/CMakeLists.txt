# -*- mode: cmake -*-

#
#  Amanzi
#   Viscous flow process kernel
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

# Define a project name
# After this command the following varaibles are defined
#   MECHANICS_SOURCE_DIR
#   MECHANICS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MECHANICS)

# External (TPL) include directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})
# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${DATA_STRUCTURES_SOURCE_DIR})
include_directories(${FLOW_SOURCE_DIR})
include_directories(${FUNCS_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})
include_directories(${MESH_SOURCE_DIR})
include_directories(${MFUNCS_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${OPERATORS_SOURCE_DIR})
include_directories(${PKS_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${STATE_SOURCE_DIR})
include_directories(${TIME_INTEGRATION_SOURCE_DIR})
include_directories(${WHETSTONE_SOURCE_DIR})

#
# Library: mechanics
#
set(mechanics_src_files Mechanics_PK.cc
                        MechanicsElasticity_PK.cc MechanicsElasticity_TI.cc
                        MechanicsSmallStrain_PK.cc MechanicsSmallStrain_TI.cc
                        MechanicsBoundaryFunction.cc
                        HydrostaticStressEvaluator.cc
                        VolumetricStrainEvaluator.cc
                        ShearModulusEvaluator.cc
   )
set(mechanics_tpl_libs ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${SEACAS_LIBRARIES})
if (Ifpack_ENABLE_HYPRE)
  list(APPEND mechanics_tpl_libs ${Ifpack_LIBRARIES} ${HYPRE_LIBRARY})
endif()
list(REVERSE mechanics_tpl_libs)
list(REMOVE_DUPLICATES mechanics_tpl_libs)
list(REVERSE mechanics_tpl_libs)

add_amanzi_library(mechanics SOURCE ${mechanics_src_files}
                   LINK_LIBS functions mesh_functions time_integration data_structures state
                             pks flow operators ${mechanics_tpl_libs})

#
# Install Targets
#
set(mechanics_inc_files 
  MechanicsElasticity_PK.hh
  MechanicsSmallStrain_PK.hh
  MechanicsBoundaryFunction.hh
  HydrostaticStressEvaluator.hh
  VolumetricStrainEvaluator.hh
)

add_install_include_file(${mechanics_inc_files})


if (BUILD_TESTS) 
  # Add UnitTest include directory
  include_directories(${UnitTest_INCLUDE_DIRS})
  include_directories(${MESH_FACTORY_SOURCE_DIR})

  # Copy test directory files if any out of source build
  if (NOT (${MECHANICS_SOURCE_DIR} EQUAL ${MECHANICS_BINARY_DIR}) )
    # execute_process(COMMAND ${CMAKE_COMMAND} -E 
    #   copy_directory ${MECHANICS_SOURCE_DIR}/test ${MECHANICS_BINARY_DIR}/test) 
    file(GLOB DataFiles "${MECHANICS_SOURCE_DIR}/test/*.xml"
                        "${MECHANICS_SOURCE_DIR}/test/*.exo")
    file(COPY ${DataFiles} DESTINATION ${MECHANICS_BINARY_DIR}/test/)
  endif()

  # Add the flow directory to the include paths
  include_directories(${MECHANICS_SOURCE_DIR})

  set(amanzi_libs geometry mesh mesh_factory state mechanics pks operators whetstone data_structures)

  # Test: 2D elastic deformation  
  add_amanzi_test(mechanics_elasticity_2D mechanics_elasticity_2D
                  KIND int
                  SOURCE test/Main.cc test/mechanics_elasticity_2D.cc
                  LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

  # Test: symmetry preservation
  add_amanzi_test(mechanics_symmetry mechanics_symmetry
                  KIND int
                  SOURCE test/Main.cc test/mechanics_symmetry.cc
                  LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

  # Test: 3D clamped beam
  add_amanzi_test(mechanics_clamped_beam mechanics_clamped_beam
                  KIND int
                  NPROCS 2
                  SOURCE test/Main.cc test/mechanics_clamped_beam.cc
                  LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

  # Test: hydrostatic stress
  add_amanzi_test(mechanics_hydrostatic_stress mechanics_hydrostatic_stress
                  KIND int
                  NPROCS 2
                  SOURCE test/Main.cc test/mechanics_hydrostatic_stress.cc
                  LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})

  # Test: spring test
  add_amanzi_test(mechanics_springs mechanics_springs
                  KIND int
                  SOURCE test/Main.cc test/mechanics_springs.cc
                  LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${HDF5_LIBRARIES})
endif()

