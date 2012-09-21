# -*- mode: cmake -*-

#
#  Amanzi
#   DivGrad operators
#
# MFD matrix methods for div-grad elliptic operators

include_directories(${ATS_SOURCE_DIR}/src/)
include_directories(upwind_scheme)
include_directories(upwind_scheme)

add_library(divgrad matrix_mfd.cc
                    upwind_scheme/upwind_cell_centered.cc
                    upwind_scheme/upwind_arithmetic_mean.cc
                    upwind_scheme/upwind_total_flux.cc
                    upwind_scheme/upwind_potential_difference.cc
                    upwind_scheme/upwind_gravity_flux.cc)
