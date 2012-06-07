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

set(PETSC_DIR /home/marc/src/ASCEM/amanzi/build/Release/TPLs/petsc-3.2)
set(PETSC_ARCH arch-linux2-c-debug)
set(PETSC_LIB_DIR ${PETSC_DIR}/${PETSC_ARCH}/lib)
#set(PETSC_LIBS petsc petscts petscsnes petscksp petscdm petscmat petscvec petscsys petscksp_single petscmat_single petscvec_single  petscsys_single)
set(PETSC_LIBS petsc)
set(PETSC_EXT_LIB_DIR /usr/local/lib /usr/lib/gcc/x86_64-linux-gnu/4.4.3 /usr/lib/x86_64-linux-gnu)

link_directories(${PETSC_EXT_LIB_DIR})
link_directories(${PETSC_LIB_DIR})
set (PETSC_EXT_LIBS X11 pthread lapack blas m dl mpich opa mpl rt pthread gcc_s mpichf90 gfortran m dl mpich opa mpl rt pthread gcc_s dl)
link_directories(${PETSC_LIB_DIR})
set(ssc_structured_grid_libs ${ssc_structured_grid_libs} ${PETSC_LIBS} ${PETSC_EXT_LIBS})


else()
  set(ssc_structured_grid_libs)
endif()
#message(">>>>> ssc_structured_grid_libs = ${ssc_structured_grid_libs}")

get_property(GEOCHEM_LIBS GLOBAL PROPERTY GEOCHEM_LIBS)
set(ssc_common_libs error_handling ${GEOCHEM_LIBS} interface atk)

set(ssc_tpl_libs_incl ${ssc_unstructured_grid_libs} ${ssc_structured_grid_libs} ${ssc_common_libs})
set(ssc_tpl_libs ${amanzi_mesh_libs} ${Teuchos_LIBRARIES} ${Epetra_LIBRARIES} ${Boost_LIBRARIES})

add_executable(amanzi main.cpp)

target_link_libraries(amanzi ${ssc_tpl_libs_incl} ${ssc_tpl_libs})
set_target_properties(amanzi PROPERTIES 
                      OUTPUT_NAME amanzi
		      OUTPUT_DIRECTORY ${SSC_BINARY_DIR})
#message(">>>>>> ssc_tpl_libs_incl = ${ssc_tpl_libs_incl}")
add_install_binary(amanzi)


if (BUILD_TESTS)

  if (ENABLE_Structured)

    add_amanzi_test(ssc-structured-test-0 amanzi "--xml_file=t0-structured.xml" KIND unit)
    add_amanzi_test(ssc-structured-test-1 amanzi "--xml_file=t1-structured.xml" KIND unit)
    add_amanzi_test(ssc-structured-test-2 amanzi "--xml_file=t2-structured.xml" KIND unit)

    if (NOT (CMAKE_SOURCE_DIR EQUAL CMAKE_BINARY_DIR))
      execute_process(COMMAND ${CMAKE_COMMAND} -E
                      copy ${SSC_SOURCE_DIR}/tests_unit/probin ${SSC_BINARY_DIR})

      execute_process(COMMAND ${CMAKE_COMMAND} -E
                      copy ${SSC_SOURCE_DIR}/tests_unit/Color.txt ${SSC_BINARY_DIR})

      execute_process(COMMAND ${CMAKE_COMMAND} -E
                      copy ${SSC_SOURCE_DIR}/tests_unit/t0-structured.xml ${SSC_BINARY_DIR})
      execute_process(COMMAND ${CMAKE_COMMAND} -E
                      copy ${SSC_SOURCE_DIR}/tests_unit/t1-structured.xml ${SSC_BINARY_DIR})
      execute_process(COMMAND ${CMAKE_COMMAND} -E
                      copy ${SSC_SOURCE_DIR}/tests_unit/t2-structured.xml ${SSC_BINARY_DIR})

    endif()
  endif()

endif()
