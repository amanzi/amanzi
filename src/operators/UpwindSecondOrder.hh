/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#ifndef AMANZI_UPWIND_SECOND_ORDER_HH_
#define AMANZI_UPWIND_SECOND_ORDER_HH_

#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "Mesh.hh"
#include "WhetStoneMeshUtils.hh"

// Operators
#include "Upwind.hh"

namespace Amanzi {
namespace Operators {

class UpwindSecondOrder : public Upwind {
 public:
  UpwindSecondOrder(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Upwind(mesh) {};
  ~UpwindSecondOrder() {};

  // main methods
  // -- initialization of control parameters
  void Init(Teuchos::ParameterList& plist);

  // -- upwind of a given cell-centered field on mesh faces
  // -- not all input parameters are use by some algorithms
  void Compute(const CompositeVector& flux,
               const std::vector<int>& bc_model,
               CompositeVector& field);

  // -- returns combined map for the original and upwinded fields.
  // -- Currently, composite vector cannot be extended on a fly.
  Teuchos::RCP<CompositeVectorSpace> Map()
  {
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1)
      ->AddComponent("grad", AmanziMesh::Entity_kind::CELL, mesh_->getSpaceDimension());
    return cvs;
  }

 private:
  int method_, order_;
  double tolerance_;
};

} // namespace Operators
} // namespace Amanzi

#endif
