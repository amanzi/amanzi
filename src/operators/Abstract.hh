/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  An abstract operator uses factory of mimetic schemes and standard
  interface for creating stiffness, mass and divergence matrices.
*/

#ifndef AMANZI_OPERATOR_ABSTRACT_HH_
#define AMANZI_OPERATOR_ABSTRACT_HH_

#include <string>
#include <vector>

// Amanzi
#include "Polynomial.hh"
#include "Tensor.hh"
#include "WhetStone_typedefs.hh"

// Operators
#include "PDE_Helper.hh"

namespace Amanzi {
namespace Operators {

class Abstract : public PDE_Helper {
 public:
  Abstract(Teuchos::ParameterList& plist, Teuchos::RCP<Operator> global_op) :
      PDE_Helper(global_op) {
    Init_(plist);
  }

  Abstract(Teuchos::ParameterList& plist, Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      PDE_Helper(mesh) {
    global_op_ = Teuchos::null;
    Init_(plist);
  }

  // main members 
  // -- required by the interface
  using PDE_Helper::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p);
  
  // -- setup
  void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K) { K_ = K; }
  void Setup(const Teuchos::RCP<const std::vector<WhetStone::Polynomial> >& Kpoly) { Kpoly_ = Kpoly; }
  void Setup(const Teuchos::RCP<const std::vector<WhetStone::VectorPolynomial> >& Kvec) { Kvec_ = Kvec; }

  // optional calculation of flux from potential p
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override {};

 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;
  Teuchos::RCP<const std::vector<WhetStone::Polynomial> > Kpoly_;
  Teuchos::RCP<const std::vector<WhetStone::VectorPolynomial> > Kvec_;

  Schema global_schema_row_, global_schema_col_;
  Schema local_schema_col_, local_schema_row_;

 private:
  void Init_(Teuchos::ParameterList& plist);

 private:
  std::string method_, matrix_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

