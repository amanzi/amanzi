# -*- mode: cmake -*-

project (FLOW_RICHARDS)
if (ENABLE_TESTS)
    add_subdirectory(steady-state)
    add_subdirectory(transient)
endif()
