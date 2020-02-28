# -*- mode: cmake -*-

# # Need this define. Errors from MSTK include files 
# # about MPI_COMM_WORLD. --lpritch
# add_definitions("-DMSTK_HAVE_MPI")
# add_definitions("-DDISABLE_PHYSICS")

# if (WITH_MSTK_2_20rc1_OR_NEWER)
#   add_definitions("-DMSTK_2_20rc1_OR_NEWER")
# endif ()

# if (WITH_MSTK_2_21rc1_OR_NEWER)
#   add_definitions("-DMSTK_2_21rc1_OR_NEWER")
# endif ()

#project(EXECUTABLE)

#
#  ATS
#    Executable
#
include_directories(${GEOCHEM_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})
include_directories(${MESH_LOGICAL_SOURCE_DIR})
include_directories(${MSTK_SOURCE_DIR})
include_directories(${CHEMPK_SOURCE_DIR})
include_directories(${MPC_TREE_SOURCE_DIR})
include_directories(${DBG_SOURCE_DIR})

#include_directories(${ATS_SOURCE_DIR}/constitutive_relations)
#include_directories(${ATS_SOURCE_DIR}/constitutive_relations/eos)
#include_directories(${ATS_SOURCE_DIR}/constitutive_relations/surface_subsurface_fluxes)
#include_directories(${ATS_SOURCE_DIR}/constitutive_relations/generic_evaluators)
include_directories(${ATS_SOURCE_DIR}/pks)
#include_directories(${ATS_SOURCE_DIR}/pks/mpc)
#include_directories(${ATS_SOURCE_DIR}/pks/energy)
#include_directories(${ATS_SOURCE_DIR}/pks/flow)
#include_directories(${ATS_SOURCE_DIR}/pks/deformation)
#include_directories(${ATS_SOURCE_DIR}/pks/transport)
include_directories(${ATS_SOURCE_DIR}/operators/upwinding)
include_directories(${ATS_SOURCE_DIR}/operators/advection)
include_directories(${ATS_SOURCE_DIR}/operators/deformation)

if (ENABLE_FATES)
  link_directories(/nh/u/dasvyat/Coastal/ats-fates-new/install/lib/)
  link_directories(/nh/u/dasvyat/Coastal/ParallelIO/install/lib/)
  link_directories(/nh/u/dasvyat/Coastal/pnetcdf-1.11.2/install/lib)
  set (VEG_LIB pk_veg)
endif()


include_evaluators_directories(LISTNAME REGISTER_AMANZI_STATE_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME ATS_RELATIONS_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_TRANSPORT_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_ENERGY_PKS_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_ENERGY_RELATIONS_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_FLOW_PKS_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_FLOW_RELATIONS_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_DEFORMATION_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_SURFACE_BALANCE_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_BGC_REG_INCLUDES)
include_evaluators_directories(LISTNAME ATS_MPC_REG_INCLUDES)
include_evaluators_directories(LISTNAME SED_TRANSPORT_REG_INCLUDES)

set(ats_src_files
  coordinator.cc
  ats_mesh_factory.cc
  simulation_driver.cc
  main.cc
  )

set(ats_inc_files
  coordinator.hh
  ats_mesh_factory.hh
  simulation_driver.hh
  )

set(amanzi_link_libs
  operators
  pks
  state
  whetstone
  time_integration
  solvers
  dbg       
  data_structures
  mesh
  mesh_audit
  mesh_functions
  functions
  geometry
  mesh_factory
  interface
  output
  mstk_mesh
  mesh_logical
  geochemrxns
  geochembase
  geochemutil
  geochemsolvers
  chemistry_pk
  mpc_tree
  transport
  )
  
set(ats_link_libs
  ats_operators
  ats_generic_evals
  ats_surf_subsurf
  ats_eos
  ats_pks
  ats_transport
  ats_sed_transport  
  ats_energy
  ats_energy_relations
  ats_flow
  ats_flow_relations
  ats_deform
  ats_bgc
  ats_surface_balance
  ats_mpc
  ats_mpc_relations  
  )  

# note, we can be inclusive here, because if they aren't enabled,
# these won't be defined and will result in empty strings.
set(tpl_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  ${Boost_LIBRARIES}
  ${PETSC_LIBRARIES}
  ${ALQUIMIA_LIBRARIES}
  ${PFLOTRAN_LIBRARIES}
  ${MSTK_LIBRARIES}
  ${SILO_LIBRARIES}
  ${HYPRE_LIBRARIES}
  ${HDF5_LIBRARIES}
  )

if (ENABLE_FATES)
  set(fates_link_libs -lfates -lparteh -lshared  -lpiof -lpioc -lpnetcdf -lgptl)
#  link_directories(ats /nh/u/dasvyat/Coastal/ats-fates-new/install/lib/)
endif()


set(ATS_LIBS
		     coordinator
                     mpc
                     #  pk_test_pks
                      # pk_energy_constant_temperature
                      pk_energy_advection_diffusion
		      pk_energy_two_phase
                      pk_energy_three_phase
                      pk_energy_surface_ice
                      pk_energy_base
		      pk_flow_permafrost
                      pk_flow_richards
                      pk_flow_snow_distribution
                      pk_flow_overland
                      pk_flow_overland_pressure
                      pk_flow_icy_overland
                      pk_transport
                      pk_coupled_transport
                      pk_reactive_transport
	              # pk_prescribed_deformation
#                      pk_prescribed_volumetric_deformation
                      pk_volumetric_deformation
                      pk_surface_balance_SEB
                      pk_BGC
                      pk_bases
                      #deform_relations_porosity
                      energy_relations
                      flow_relations
                      flow_relations_surface_subsurface_fluxes
                      generic_evaluators
                      relations_eos
                      relations_ewc
                      bc_factory
                      advection
                      divgrad
                      deformation_operator
)

set(AMANZI_LIBS
                      amanzi_operators
                      amanzi_pks
                      amanzi_state
	              amanzi_whetstone
                      amanzi_time_integration
	              amanzi_solvers
                      amanzi_data_structures
                      amanzi_mesh
                      amanzi_mesh_audit
                      amanzi_mesh_functions
                      amanzi_functions
                      amanzi_geometry
                      amanzi_mesh_factory
                      amanzi_interface
                      amanzi_output
                      amanzi_mstk_mesh
		      amanzi_mesh_logical
                      amanzi_geochemrxns
		      amanzi_geochembase
		      amanzi_geochemutil
		      amanzi_geochemsolvers
		      amanzi_chemistry_pk
		      amanzi_mpc_tree
		      amanzi_transport
                      amanzi_dbg
)

# set(GEOCHEM_LIBS 
#                  /home/scratch/dasvyat/Amanzi/build/amanzi_default/lib/libchemistry_pk.a 
#                  /home/scratch/dasvyat/Amanzi/build/amanzi_default/lib/libgeochembase.a 
#                  /home/scratch/dasvyat/Amanzi/build/amanzi_default/lib/libgeochemsolvers.a 
#                  /home/scratch/dasvyat/Amanzi/build/amanzi_default/lib/libgeochemutil.a 
#                  /home/scratch/dasvyat/Amanzi/build/amanzi_default/lib/libmpc_tree.a 
#                  /home/scratch/dasvyat/Amanzi/build/amanzi_default/lib/libtransport.a
# )

target_link_libraries(ats ${ATS_LIBS} ${AMANZI_LIBS})# ${GEOCHEM_LIBS})


install(TARGETS ats DESTINATION bin)

#------------------------------------------------------------------------------#
# ATS F90 test program
#------------------------------------------------------------------------------#

if (Amanzi_DISABLE_MPI_CXX_BINDINGS)  

   add_executable(ats_f90 main.f90)

   target_link_libraries(
        # Target
	ats_f90

        # ATS cross lib
	atsx

    ${ATS_LIBS}
	${AMANZI_LIBS}
	${Amanzi_TPL_Trilinos_LIBRARIES}

        # MPI
	${Amanzi_MPI_Fortran_LIBRARIES}
	${Amanzi_MPI_C_LIBRARIES}
   )

   set_target_properties(ats_f90 PROPERTIES LINKER_LANGUAGE Fortran)

   #---------------------------------------------------------------------------#
   # Install rule
   #---------------------------------------------------------------------------#

   install(TARGETS ats_f90 DESTINATION bin)



   add_executable(ats_test_clm_driver test_clm_driver.cc)

   target_link_libraries(
        # Target
	ats_test_clm_driver

        # ATS cross lib
	atsx

    ${ATS_LIBS}
	${AMANZI_LIBS}
	${Amanzi_TPL_Trilinos_LIBRARIES}

        # MPI
	${Amanzi_MPI_C_LIBRARIES}
   )

   #---------------------------------------------------------------------------#
   # Install rule
   #---------------------------------------------------------------------------#

 install(TARGETS ats_test_clm_driver DESTINATION bin)

endif()

# add_amanzi_executable(ats
#   SOURCE ${ats_src_files}
#   LINK_LIBS ${fates_link_libs} ${ats_link_libs} ${amanzi_link_libs} ${tpl_link_libs}
#   OUTPUT_NAME ats
#   OUTPUT_DIRECTORY ${ATS_BINARY_DIR})
