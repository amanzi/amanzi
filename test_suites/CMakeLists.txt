# -*- mode: cmake -*-

project(AMANZI_TESTS)

if (ENABLE_TESTS)
    add_subdirectory(verification)
endif()
