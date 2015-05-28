# -*- mode: cmake -*-

project (FLOW_RICHARDS_STEADY)

if (ENABLE_TESTS)
    add_subdirectory(infiltration_clay_sand_1d)
    add_subdirectory(infiltration_loam_sand_1d)
    add_subdirectory(infiltration_sand_loam_1d)
    add_subdirectory(unconfined_layered_2d)
    add_subdirectory(unconfined_no_recharge_1d)
    # EIB - skipped for now because it takes forever to run
    add_subdirectory(unconfined_recharge_1d)
    add_subdirectory(unconfined_seepage_1d)
endif()
