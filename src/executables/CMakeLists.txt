# -*- mode: cmake -*-

# Need this define. Errors from MSTK include files 
# about MPI_COMM_WORLD. --lpritch
add_definitions("-DMSTK_HAVE_MPI")
add_definitions("-DDISABLE_PHYSICS")

if (WITH_MSTK_2_20rc1_OR_NEWER)
  add_definitions("-DMSTK_2_20rc1_OR_NEWER")
endif ()

if (WITH_MSTK_2_21rc1_OR_NEWER)
  add_definitions("-DMSTK_2_21rc1_OR_NEWER")
endif ()

#project(EXECUTABLE)

#
#  ATS
#    Executable
#
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
include_directories(${ATS_SOURCE_DIR}/src/pks/coupled_transport)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad/upwind_scheme)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/operators/deformation)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/base)
include_directories(${ATS_SOURCE_DIR}/src/pks/transport)
include_directories(${ATS_SOURCE_DIR}/src/pks/transport/transport_amanzi)
include_directories(${ATS_SOURCE_DIR}/src/pks/chemistry)
include_directories(${ATS_SOURCE_DIR}/src/pks/chemistry/chemistry_amanzi)

include_directories(${ATS_SOURCE_DIR}/src/coordinator)

if (ENABLE_FATES)
  link_directories(/nh/u/dasvyat/Coastal/ats-fates-new/install/lib/)
  link_directories(/nh/u/dasvyat/Coastal/ParallelIO/install/lib/)
  link_directories(/nh/u/dasvyat/Coastal/pnetcdf-1.11.2/install/lib)
  set (VEG_LIB pk_veg)
endif()

include_evaluators_directories(LISTNAME ATS_POROSITY_EVAL_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_EOS_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_SURFACE_SUBSURFACE_FLUXES_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_GENERIC_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_RELATIONS_INCLUDES)
include_evaluators_directories(LISTNAME WRM_MODELS_EVALUATORS_FACTORIES_INCLUDES)
include_evaluators_directories(LISTNAME WMR_EVALUATORS_FACTORIES_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_ICY_OVERLAND_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_OVERLAND_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_OVERLAND_HEAD_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_OVERLAND_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME RICHARDS_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME PERMAFROST_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME DEFORMATION_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_ADVECTION_DIFFUSION_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_CONSTANT_TEMPERATURE_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_RELATIONS_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_SURFACE_ICE_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_TWO_PHASE_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_THREE_PHASE_REG_INCLUDES)
include_evaluators_directories(LISTNAME SURFACE_BALANCE_SEB_REG_INCLUDES)
include_evaluators_directories(LISTNAME BGC_REG_INCLUDES)
include_evaluators_directories(LISTNAME TEST_PKS_REG_INCLUDES)
include_evaluators_directories(LISTNAME MPC_PK_FAC_INCLUDES)
#include_evaluators_directories(LISTNAME TRANSPORT_PASSIVE_TRACER_REG_INCLUDES)
include_evaluators_directories(LISTNAME TRANSPORT_AMANZI_REG_INCLUDES)
  


#add_executable(ats main.cc plant_1D_mesh.cc ats_mesh_factory.cc simulation_driver.cc )
add_executable(ats main.cc ats_mesh_factory.cc simulation_driver.cc )

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
#                      pk_flow_overland
                      pk_flow_overland_pressure
                      pk_flow_icy_overland
                      pk_transport

                      pk_morphology
                      pk_coupled_transport
                      pk_sediment_transport

	              # pk_prescribed_deformation
#                      pk_prescribed_volumetric_deformation
                      pk_volumetric_deformation
                      pk_surface_balance_SEB
                      pk_BGC
                      ${VEG_LIB}
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

if (ENABLE_FATES)
  set(FATES_LIBS -lfates -lparteh -lshared  -lpiof -lpioc -lpnetcdf -lgptl)
  link_directories(ats /nh/u/dasvyat/Coastal/ats-fates-new/install/lib/)
endif()




#target_link_libraries(ats  ${ATS_LIBS} ${FATES_LIBS} ${AMANZI_LIBS} hdf5::hdf5_hl-static)
target_link_libraries(ats  ${ATS_LIBS} ${FATES_LIBS} ${AMANZI_LIBS})
install(TARGETS ats DESTINATION bin)


message("build tests = ${BUILD_TESTS}")
#if ( BUILD_TESTS )
if ( 0 )  
    include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})
    include_directories(".")

    # Test:
    # add_amanzi_test(plant_mesh plant_mesh
    #                 KIND unit
    #                 SOURCE test/main.cc test/plant_mesh.cc plant_1D_mesh.cc
    #                 LINK_LIBS ${ATS_LIBS} ${AMANZI_LIBS} ${Amanzi_TPL_UnitTest_LIBRARIES})
endif()   
