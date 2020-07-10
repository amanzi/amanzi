# -*- mode: cmake -*-

#
#  ATS
#    Lake thermo PK class
#

# additional common components
#add_subdirectory(constitutive_relations)

# ATS include directories
include_directories(${ATS_SOURCE_DIR}/pks)
include_directories(${ATS_SOURCE_DIR}/operators/advection)
#include_directories(${ATS_SOURCE_DIR}/operators/upwinding)
#include_directories(${ATS_SOURCE_DIR}/pks/energy/constitutive_relations/enthalpy)
#include_directories(${ATS_SOURCE_DIR}/pks/energy/constitutive_relations/energy)
#include_directories(${ATS_SOURCE_DIR}/pks/energy/constitutive_relations/internal_energy)
#include_directories(${ATS_SOURCE_DIR}/pks/energy/constitutive_relations/thermal_conductivity)
#include_directories(${ATS_SOURCE_DIR}/pks/energy/constitutive_relations/source_terms)
include_directories(${ATS_SOURCE_DIR}/constitutive_relations/eos)


set(ats_lake_src_files
  lake_thermo_pk.cc
  )

set(ats_lake_inc_files
  lake_thermo_pk.hh
  )


set(ats_lake_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  error_handling
  atk
  mesh
  data_structures
  whetstone
  operators
  solvers
  time_integration
  state
  pks
  ats_operators
  ats_eos
  ats_pks
 # ats_lake_relations
  )


add_amanzi_library(ats_lake
                   SOURCE ${ats_lake_src_files}
                   HEADERS ${ats_lake_inc_files}
		   LINK_LIBS ${ats_lake_link_libs})



#================================================
# register evaluators/factories/pks
register_evaluator_with_factory(
  HEADERFILE lake_thermo_pk_reg.hh
  LISTNAME   ATS_LAKE_PKS_REG
  )

#generate_evaluators_registration_header(
#  HEADERFILE lake_thermo_pk_reg.hh
#  LISTNAME   ATS_LAKE_PKS_REG
#  INSTALL    True
#  )

