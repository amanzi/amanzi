# -*- mode: cmake -*-

#
#  Amanzi
#   DivGrad operators
#
# MFD matrix methods for div-grad elliptic operators

add_library(deformation_operator MatrixVolumetricDeformation.cc
                                 Matrix_PreconditionerDelegate.cc)

install(TARGETS deformation_operator DESTINATION lib)
