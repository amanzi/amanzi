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

/*
This is the discretization package, release beta.
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "tensor.hh"
#include "mfd3d.hh"


namespace Amanzi {
namespace WhetStone {

class NLFV { 
 public:
  NLFV(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh) {};
  ~NLFV() {};

  void HarmonicAveragingPoint(int face, std::vector<Tensor>& T,
                              AmanziGeometry::Point& hap, double& hap_weight);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

