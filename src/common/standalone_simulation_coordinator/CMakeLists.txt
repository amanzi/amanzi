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
include(BinaryManager)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})

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
  include_directories(${STRUCTUREDUTILS_SOURCE_DIR})
  link_directories(${STRUCTUREDUTILS_SOURCE_DIR})
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
  include_directories(${MESH_FACTORY_SOURCE_DIR})
  include_directories(${GEOMETRY_SOURCE_DIR})
  include_directories(${OUTPUT_SOURCE_DIR})
  include_directories(${STATE_SOURCE_DIR})
  include_directories(${TRANSPORT_SOURCE_DIR})
  include_directories(${FLOW_SOURCE_DIR})
  include_directories(${MPC_SOURCE_DIR})
  include_directories(${EXO_SOURCE_DIR})
  include_directories(${EXEC_SOURCE_DIR})

  link_directories(${EXEC_SOURCE_DIR})

# The mesh frameworks
  set(amanzi_mesh_libs "mesh_factory")

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
  set(ssc_unstructured_grid_libs udriv mpc state transport flow chemistry_pk whetstone matrix_audit)
else()
  set(ssc_unstructured_grid_libs)
endif()

if(ENABLE_Structured)
  set(structured_grid_libs pmamr gsint boxprob boxdriv structutils ${CCSE_LIBRARIES})
  link_directories(${CCSE_EXT_LIBRARY_DIRS})
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
  list(APPEND ssc_structured_grid_libs ${PETSC_LIBRARIES})
else()
  set(ssc_structured_grid_libs)
endif()
#message(">>>>> ssc_structured_grid_libs = ${ssc_structured_grid_libs}")

get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)
set(ssc_common_libs error_handling ${GEOCHEM_LIBS} interface atk)

set(ssc_tpl_libs_incl ${ssc_unstructured_grid_libs} ${ssc_structured_grid_libs} ${ssc_common_libs})
set(ssc_tpl_libs ${amanzi_mesh_libs} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${Boost_LIBRARIES})

if (ENABLE_Structured)
  list(APPEND ssc_tpl_libs ${CCSE_EXT_LIBRARIES})
endif()


# 
# Main: Amanzi
#
#print_variable(ssc_tpl_libs_incl)
#print_variable(ssc_tpl_libs)
#message(">>>>>> ssc_tpl_libs_incl = ${ssc_tpl_libs_incl}")
add_amanzi_executable(amanzi
                      SOURCE main.cpp
                      LINK_LIBS ${ssc_tpl_libs_incl} ${ssc_tpl_libs}
                      OUTPUT_NAME amanzi
                      OUTPUT_DIRECTORY ${SSC_BINARY_DIR})

if (BUILD_TESTS AND ENABLE_Structured) 

    # Copy test directory files if an out of source build
    if (NOT (${SSC_SOURCE_DIR} EQUAL ${SSC_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${SSC_SOURCE_DIR}/test ${SSC_BINARY_DIR}/test) 
    endif()

    add_amanzi_test(ssc-structured-transport-0 amanzi "--xml_file=test/t0-structured.xml" KIND unit)
    add_amanzi_test(ssc-structured-flow-1 amanzi "--xml_file=test/t1-structured.xml" KIND unit)

endif()

