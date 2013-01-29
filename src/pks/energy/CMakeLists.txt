# -*- mode: cmake -*-

#
#  ATS
#    Flow PK class
#

# ATS include directories

# additional common components
add_subdirectory(bc_factory)
add_subdirectory(constitutive_relations)

# PKs
add_subdirectory(constant_temperature)
add_subdirectory(advection_diffusion)
add_subdirectory(base)
add_subdirectory(two_phase)
add_subdirectory(three_phase)


