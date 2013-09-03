/*
This is the mimetic discretization component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Version: 2.0
Release name: naka-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#ifndef __NLFV_HH__
#define __NLFV_HH__

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "WhetStone_typedefs.hh"
#include "tensor.hh"


namespace Amanzi {
namespace WhetStone {

class NLFV { 
 public:
  NLFV(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh) {};
  ~NLFV() {};

  void HarmonicAveragingPoint(int face, std::vector<Tensor>& T,
                              AmanziGeometry::Point& p, double& weight);

  void MaximumDecomposition(const AmanziGeometry::Point& conormal, 
                            const std::vector<AmanziGeometry::Point>& tau,
                            double* w1, double* w2, int* i1, int* i2);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

