/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Nonlinear finite volume method.
*/

#ifndef AMANZI_WHETSTONE_NLFV_HH_
#define AMANZI_WHETSTONE_NLFV_HH_

#include "Teuchos_RCP.hpp"

#include "MeshLight.hh"
#include "Point.hh"

#include "WhetStoneDefs.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class NLFV {
 public:
  NLFV() : mesh_(Teuchos::null){};
  NLFV(Teuchos::RCP<const AmanziMesh::MeshLight> mesh) : mesh_(mesh){};
  ~NLFV(){};

  void HarmonicAveragingPoint(int f,
                              int c1,
                              int c2,
                              const AmanziGeometry::Point& Tn1,
                              const AmanziGeometry::Point& Tn2,
                              AmanziGeometry::Point& p,
                              double& weight);

  int PositiveDecomposition(int id1,
                            const std::vector<AmanziGeometry::Point>& tau,
                            const AmanziGeometry::Point& conormal,
                            double* ws,
                            int* ids);

 private:
  Teuchos::RCP<const AmanziMesh::MeshLight> mesh_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
