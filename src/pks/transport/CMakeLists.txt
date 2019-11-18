# -*- mode: cmake -*-

#
#  ATS
#    Transport Amanzi
#


# ATS include directories
add_subdirectory(sediment_transport)

#================================================
# register evaluators/factories/pks

set(ats_transport_inc_files
  TransportDefs.hh
  MDM.hh
  MDM_Isotropic.hh
  MDMFactory.hh
  MDM_LichtnerKelkarRobinson.hh
  MDMPartition.hh
  MDM_Bear.hh
  MDM_BurnettFrind.hh
  DiffusionPhase.hh
  MaterialProperties.hh
  MultiscaleTransportPorosityPartition.hh
  TransportDomainFunction.hh  
  MultiscaleTransportPorosity_DPM.hh
  TransportSourceFunction.hh
  TransportSourceFunction_Alquimia.hh
  MultiscaleTransportPorosity_GDPM.hh
  MultiscaleTransportPorosity.hh
  TransportBoundaryFunction.hh
  MultiscaleTransportPorosityFactory.hh
  TransportBoundaryFunction_Alquimia.hh
  Transport_PK_ATS.hh
  )


set(ats_transport_link_libs
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
  chempk
  ats_operators
  ats_eos
  ats_pks
  )

register_evaluator_with_factory(
  HEADERFILE sediment_transport/sediment_transport_reg.hh
  LISTNAME   TRANSPORT_AMANZI_REG
  )

add_amanzi_library(ats_transport
                   SOURCE ${ats_transport_src_files}
                   HEADERS ${ats_transport_inc_files}
		   LINK_LIBS ${ats_transport_link_libs})

#sediment diposition and erosion rates
register_evaluator_with_factory(
  HEADERFILE sediment_transport/erosion_evaluator_reg.hh
  LISTNAME SEDIMENT_TRANSPORT_REGISTRATIONS
  )
register_evaluator_with_factory(
  HEADERFILE sediment_transport/settlement_evaluator_reg.hh
  LISTNAME SEDIMENT_TRANSPORT_REGISTRATIONS
  )
register_evaluator_with_factory(
  HEADERFILE sediment_transport/trapping_evaluator_reg.hh
  LISTNAME SEDIMENT_TRANSPORT_REGISTRATIONS
  )
register_evaluator_with_factory(
  HEADERFILE sediment_transport/organic_matter_evaluator_reg.hh
  LISTNAME SEDIMENT_TRANSPORT_REGISTRATIONS
  )

generate_evaluators_registration_header( 
  HEADERFILE sediment_transport_registration.hh   
  LISTNAME   SEDIMENT_TRANSPORT_REGISTRATIONS
  INSTALL    True
  )

#================================================
# register evaluators/factories/pks

# transport_amanzi
register_evaluator_with_factory(
  HEADERFILE Transport_PK_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )

# -- mechamical dispersion models registration
register_evaluator_with_factory(
  HEADERFILE MDMFactory_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE MDM_Isotropic_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE MDM_Bear_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE MDM_BurnettFrind_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE MDM_LichtnerKelkarRobinson_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )

# -- multiscale models registration
register_evaluator_with_factory(
  HEADERFILE MultiscaleTransportPorosityFactory_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE MultiscaleTransportPorosity_DPM_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )
register_evaluator_with_factory(
  HEADERFILE MultiscaleTransportPorosity_GDPM_reg.hh
  LISTNAME   ATS_TRANSPORT_REG
  )

generate_evaluators_registration_header(
  HEADERFILE ats_transport_registration.hh
  LISTNAME   ATS_TRANSPORT_REG
  INSTALL    True
  )

