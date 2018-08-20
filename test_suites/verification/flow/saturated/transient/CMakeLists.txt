# -*- mode: cmake -*-

project (FLOW_SATURATED_TRANSIENT)

if (ENABLE_TESTS)
    add_subdirectory(boundedDomain_2d)
    add_subdirectory(butler_pod_2d)
    add_subdirectory(butler_strip_2d)
    add_subdirectory(hantush_anisotropic_2d)
    add_subdirectory(theis_isotropic_1d)
endif()
