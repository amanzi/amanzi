/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

Elasticity operator on meshes with curved faces.

*/

#ifndef AMANZI_OPERATOR_PDE_ELASTICITY_CURVED_FACE_HH_
#define AMANZI_OPERATOR_PDE_ELASTICITY_CURVED_FACE_HH_

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "exceptions.hh"
#include "CompositeVector.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "BilinearForm.hh"
#include "MFD3D_ElasticityWeakSym_CurvedFace.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_DiffusionCurvedFace.hh"
#include "PDE_HelperDiscretization.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class PDE_ElasticityCurvedFace : public PDE_Abstract {
 public:
  PDE_ElasticityCurvedFace(Teuchos::ParameterList& plist,
                           const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_Abstract(plist, mesh)
  {
    Init_(plist);
  };

 private:
  void Init_(Teuchos::ParameterList& plist)
  {
    PDE_DiffusionCurvedFace pde(plist, mesh_);
    bf_ = pde.get_bf();
    Teuchos::rcp_dynamic_cast<WhetStone::MFD3D_ElasticityWeakSym_CurvedFace>(mfd_)->set_generalized_centroids(bf_);
  }

 protected:
  std::shared_ptr<std::vector<AmanziGeometry::Point>> bf_;
};

} // namespace Operators
} // namespace Amanzi

#endif
