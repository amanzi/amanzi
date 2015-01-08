# -*- mode: cmake -*-

project (VERIFICATION)

if (ENABLE_TESTS)
    add_subdirectory(flow)
    add_subdirectory(transport)
endif()
