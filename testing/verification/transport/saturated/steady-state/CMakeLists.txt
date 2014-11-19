# -*- mode: cmake -*-

project (TRANSPORT_SATURATED_STEADY)

if (ENABLE_TESTS)
    add_subdirectory(dispersion_45_point_2d)
    add_subdirectory(dispersion_aligned_point_2d)
endif()
