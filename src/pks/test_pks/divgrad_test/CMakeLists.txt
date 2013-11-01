# -*- mode: cmake -*-
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/factory)
include_directories(${ATS_SOURCE_DIR}/src/pks/bc_factory)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)

add_library(pk_test_divgrad divgrad_test.cc)

