# -*- mode: cmake -*-

project (FLOW_TESTS)

if (ENABLE_TESTS)
    add_subdirectory(richards)
    add_subdirectory(saturated)
endif()
