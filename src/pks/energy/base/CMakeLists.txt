# -*- mode: cmake -*-
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad/upwind_scheme)

add_library(pk_energy_base energy_base_pk.cc energy_base_ti.cc energy_base_physics.cc)

install(TARGETS pk_energy_base DESTINATION lib)


