/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Examples of usage this operator are in test/operators_elasticity.cc
  and test/operators_stokes.cc
*/

#ifndef AMANZI_OPERATOR_PDE_ELASTICITY_HH_
#define AMANZI_OPERATOR_PDE_ELASTICITY_HH_

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "exceptions.hh"
#include "Tensor.hh"
#include "CompositeVector.hh"

// Amanzi::Operators
#include "BilinearForm.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_HelperDiscretization.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class PDE_Elasticity : public PDE_HelperDiscretization {
 public:
  PDE_Elasticity(Teuchos::ParameterList& plist, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_HelperDiscretization(mesh), K_(Teuchos::null), K_default_(1.0)
  {
    global_op_ = Teuchos::null;
    pde_type_ = PDE_ELASTICITY;
    Init_(plist);
  }

  // main virtual members
  // -- setup
  void SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor>>& K);
  void SetTensorCoefficient(double K);

  // -- creation of an operator
  using PDE_HelperDiscretization::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p) override;

  // -- postprocessing: calculated stress u from displacement p
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override{};

 protected:
  void Init_(Teuchos::ParameterList& plist);

 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K_;
  double K_default_;

  Teuchos::RCP<WhetStone::BilinearForm> mfd_;
  AmanziMesh::Entity_kind base_;
};

} // namespace Operators
} // namespace Amanzi

#endif
