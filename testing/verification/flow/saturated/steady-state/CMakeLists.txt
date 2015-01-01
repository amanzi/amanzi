# -*- mode: cmake -*-

project (FLOW_SATURATED_STEADY)

if (ENABLE_TESTS)
    add_subdirectory(linear_flux_head_1d)
    add_subdirectory(linear_head_head_1d)
    add_subdirectory(linear_materials_parallel_1d)
    add_subdirectory(linear_materials_serial_1d)
endif()
