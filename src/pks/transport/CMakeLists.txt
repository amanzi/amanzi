# -*- mode: cmake -*-

#
#  ATS
#    Transport Amanzi
#

# ATS include directories
#add_subdirectory(passive_tracer)
add_subdirectory(transport_amanzi)

#================================================
# register evaluators/factories/pks

# passive_tracer

register_evaluator_with_factory(
  HEADERFILE passive_tracer/passive_tracer_reg.hh
  LISTNAME   TRANSPORT_PASSIVE_TRACER_REG
  )

generate_evaluators_registration_header(
  HEADERFILE transport_passive_tracer_registration.hh
  LISTNAME   TRANSPORT_PASSIVE_TRACER_REG
  INSTALL    True
  )

# transport_amanzi

register_evaluator_with_factory(
  HEADERFILE transport_amanzi/Transport_PK_reg.hh
  LISTNAME   TRANSPORT_AMANZI_REG
  )


generate_evaluators_registration_header(
  HEADERFILE transport_amanzi_registration.hh
  LISTNAME   TRANSPORT_AMANZI_REG
  INSTALL    True
  )


# -- mechamical dispersion models registration
register_evaluator_with_factory(
  HEADERFILE transport_amanzi/MDMFactory_reg.hh
  LISTNAME   MDM_TRANSPORT_REGISTRATIONS
  )

register_evaluator_with_factory(
  HEADERFILE transport_amanzi/MDM_Isotropic_reg.hh
  LISTNAME   MDM_TRANSPORT_REGISTRATIONS
  )

register_evaluator_with_factory(
  HEADERFILE transport_amanzi/MDM_Bear_reg.hh
  LISTNAME   MDM_TRANSPORT_REGISTRATIONS
  )

register_evaluator_with_factory(
  HEADERFILE transport_amanzi/MDM_BurnettFrind_reg.hh
  LISTNAME   MDM_TRANSPORT_REGISTRATIONS
  )

register_evaluator_with_factory(
  HEADERFILE transport_amanzi/MDM_LichtnerKelkarRobinson_reg.hh
  LISTNAME   MDM_TRANSPORT_REGISTRATIONS
  )

generate_evaluators_registration_header( 
  HEADERFILE mdm_transport_registration.hh   
  LISTNAME   MDM_TRANSPORT_REGISTRATIONS
  INSTALL    True
  )

# -- multiscale models registration
register_evaluator_with_factory(
  HEADERFILE transport_amanzi/MultiscaleTransportPorosityFactory_reg.hh
  LISTNAME   MSM_TRANSPORT_REGISTRATIONS
  )

register_evaluator_with_factory(
  HEADERFILE transport_amanzi/MultiscaleTransportPorosity_DPM_reg.hh
  LISTNAME   MSM_TRANSPORT_REGISTRATIONS
  )

generate_evaluators_registration_header( 
  HEADERFILE multiscale_transport_registration.hh   
  LISTNAME   MSM_TRANSPORT_REGISTRATIONS
  INSTALL    True
  )

