# -*- mode: cmake -*-
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/surface_subsurface_fluxes)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/overland_conductivity)
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/elevation)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad/upwind_scheme)
include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})

add_library(pk_flow_overland_pressure
        overland_pressure_pk.cc overland_pressure_ti.cc overland_pressure_physics.cc
        height_evaluator.cc
	volumetric_height_evaluator.cc	
        effective_height_evaluator.cc
        overland_pressure_water_content_evaluator.cc
	overland_subgrid_water_content_evaluator.cc)

install(TARGETS pk_flow_overland_pressure DESTINATION lib)
