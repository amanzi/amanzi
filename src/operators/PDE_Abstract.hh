/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_PDE_ABSTRACT_HH_
#define AMANZI_OPERATOR_PDE_ABSTRACT_HH_

#include <string>
#include <vector>

// Amanzi
#include "BilinearForm.hh"
#include "Polynomial.hh"
#include "Tensor.hh"
#include "VectorPolynomial.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "PDE_HelperDiscretization.hh"

namespace Amanzi {
namespace Operators {

class PDE_Abstract : public PDE_HelperDiscretization {
 public:
  PDE_Abstract(Teuchos::ParameterList& plist, Teuchos::RCP<Operator> global_op)
    : PDE_HelperDiscretization(global_op)
  {
    Init_(plist);
  }

  PDE_Abstract(Teuchos::ParameterList& plist,
               Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : PDE_HelperDiscretization(mesh)
  {
    global_op_ = Teuchos::null;
    Init_(plist);
  }
  ~PDE_Abstract(){};

  // main members
  // -- required by the interface
  using PDE_HelperDiscretization::UpdateMatrices;
  virtual void
  UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                 const Teuchos::Ptr<const CompositeVector>& p) override;

  // -- setup
  void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor>>& K) { K_ = K; }
  void
  SetupPoly(const Teuchos::RCP<const std::vector<WhetStone::Polynomial>>& Kpoly)
  {
    Kpoly_ = Kpoly;
  }
  void SetupPolyVector(
    const Teuchos::RCP<const std::vector<WhetStone::VectorPolynomial>>& Kvec)
  {
    Kvec_ = Kvec;
  }

  // optional calculation of flux from potential p
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override{};

 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K_;
  Teuchos::RCP<const std::vector<WhetStone::Polynomial>> Kpoly_;
  Teuchos::RCP<const std::vector<WhetStone::VectorPolynomial>> Kvec_;

  Schema global_schema_row_, global_schema_col_;
  Schema local_schema_col_, local_schema_row_;

 private:
  void Init_(Teuchos::ParameterList& plist);

 private:
  std::string matrix_;
  bool grad_on_test_;

  Teuchos::RCP<WhetStone::BilinearForm> mfd_;
};

} // namespace Operators
} // namespace Amanzi

#endif
