# -*- mode: cmake -*-
#
# Define a project name
# After this command the following varaibles are defined
#   SSC_SOURCE_DIR
#   SSC_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(SSC)

include(CCSEOptions)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})

# The interfaces
include_directories(${INTERFACE_PLATFORM_DIR})

#
# The structured grid pieces
#
if ( ENABLE_Structured )
  include_directories(${CCSE_INCLUDE_DIRS})
  include_directories(${PMAMR_INCLUDE_DIRS})
  link_directories(${CCSE_LIBRARY_DIR})
  include_directories(${BOXDRIV_SOURCE_DIR})
  link_directories(${BOXDRIV_SOURCE_DIR})
  add_definitions("-DENABLE_Structured")
endif()

#
# The common pieces
#
include_directories(${DBC_SOURCE_DIR})
include_directories(${CHEM_INCLUDES_DIR})

#
# The unstructured grid pieces
#
if ( ENABLE_Unstructured )
  include_directories(${ATK_SOURCE_DIR})
  include_directories(${MESH_SOURCE_DIR})
  include_directories(${MESH_FACTORY_SOURCE_DIR})
  include_directories(${GEOMETRY_SOURCE_DIR})
  include_directories(${OUTPUT_SOURCE_DIR})
  include_directories(${STATE_SOURCE_DIR})
  include_directories(${TRANSPORT_SOURCE_DIR})
  include_directories(${FLOW_SOURCE_DIR})
  include_directories(${MPC_SOURCE_DIR})
  include_directories(${MESH_SIMPLE_SOURCE_DIR})
  include_directories(${MESH_DATA_SOURCE_DIR})
  include_directories(${EXO_SOURCE_DIR})
  include_directories(${EXEC_SOURCE_DIR})

  link_directories(${EXEC_SOURCE_DIR})

# The mesh frameworks
  set(amanzi_mesh_libs "")
  if ( ENABLE_STK_Mesh )
      include_directories(${STK_INCLUDE_DIRS})
      include_directories(${STK_SOURCE_DIR})
      list(APPEND amanzi_mesh_libs stk_mesh exodus)
      add_definitions("-DHAVE_STK_MESH")
  endif()
  if(ENABLE_MOAB_Mesh)
      include_directories(${MOAB_INCLUDE_DIRS})
      include_directories(${MOAB_SOURCE_DIR})
      list(APPEND amanzi_mesh_libs moab_mesh)
      add_definitions("-DHAVE_MOAB_MESH")
  endif()    
  if(ENABLE_MSTK_Mesh)
      include_directories(${MSTK_INCLUDE_DIRS})
      include_directories(${MSTK_SOURCE_DIR})
      list(APPEND amanzi_mesh_libs mstk_mesh)
  endif()    

  add_definitions("-DENABLE_Unstructured")
endif()

# External (TPL) include directories
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})

#
# More unstructured grid pieces
#
if ( ENABLE_Unstructured )
  if(ENABLE_CGNS)
      include_directories(${CGNS_INCLUDE_DIRS})
  endif()    

  # Need to move this up to the root -- lpritch
  add_definitions("-DUSE_MPI")
endif()


#
# Binary: amanzi-demo-driver
#
if ( ENABLE_Unstructured )
  set(ssc_unstructured_grid_libs udriv mpc state transport flow chemistry)
else()
  set(ssc_unstructured_grid_libs)
endif()

if(ENABLE_Structured)
  set(structured_grid_libs pmamr gsint boxprob boxdriv ${CCSE_LIBRARIES})
  if (${CMAKE_SYSTEM_NAME} STREQUAL Darwin)
    #message(">>>>>> using darwin ld flags")
    set(ld_flag_prefix -force_load)
    set(ld_flag_suffix)
  elseif (${CMAKE_SYSTEM_NAME} STREQUAL Linux)
    #message(">>>>>> using linux ld flags")
    set(ld_flag_prefix -Wl,-whole-archive)
    set(ld_flag_suffix -Wl,-no-whole-archive)
  else()
    message(FATAL_ERROR "Linking of structured grid is not supported on this system: '${CMAKE_SYSTEM_NAME}'")
  endif()
  foreach(lib ${structured_grid_libs})
    set(ssc_structured_grid_libs ${ssc_structured_grid_libs} ${ld_flag_prefix} ${lib} ${ld_flag_suffix})
  endforeach(lib)
else()
  set(ssc_structured_grid_libs)
endif()
#message(">>>>> ssc_structured_grid_libs = ${ssc_structured_grid_libs}")

get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)
set(ssc_common_libs error_handling ${GEOCHEM_LIBS})

set(ssc_tpl_libs_incl ${ssc_unstructured_grid_libs} ${ssc_structured_grid_libs} ${ssc_common_libs})
set(ssc_tpl_libs ${amanzi_mesh_libs} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${Boost_LIBRARIES})

add_executable(ssc main.cpp)

target_link_libraries(ssc ${ssc_tpl_libs_incl} ${ssc_tpl_libs})
#message(">>>>>> ssc_tpl_libs_incl = ${ssc_tpl_libs_incl}")

if (BUILD_TESTS)
endif()





