/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discontinuous Galerkin method.
*/

#ifndef AMANZI_WHETSTONE_DG_HH_
#define AMANZI_WHETSTONE_DG_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "WhetStone_typedefs.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class DG { 
 public:
  DG() : mesh_(Teuchos::null) {};
  DG(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh) {};
  ~DG() {};

  int TaylorMassMatrix(int c);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

