# -*- mode: cmake -*-

project (FLOW_RICHARDS_TRANSIENT)

if (ENABLE_TESTS)
    add_subdirectory(infiltration_dry_sand_1d)
endif()
