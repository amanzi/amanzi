/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  An abstract operator uses factory of mimetic schemes and standard
  interface for creating stiffness, mass and divergence matrices.

  Examples of usage this operator are in test/operators_stokes.cc
  and test/operators_diffusion_curved.cc
  In the first example, we set up a discrete divergence operator
  that corersponds to a rectangular matrix. In the second example,
  we set up an elliptic operator when Hermite-type degrees of 
  freedom are used on curved faces.
*/

#ifndef AMANZI_OPERATOR_PDE_ABSTRACT_HH_
#define AMANZI_OPERATOR_PDE_ABSTRACT_HH_

#include <string>
#include <vector>

// Amanzi
#include "BilinearForm.hh"
#include "Polynomial.hh"
#include "Tensor.hh"
#include "VectorObjects.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "PDE_HelperDiscretization.hh"

namespace Amanzi {
namespace Operators {

class PDE_Abstract : public PDE_HelperDiscretization {
 public:
  PDE_Abstract(Teuchos::ParameterList& plist, Teuchos::RCP<Operator> global_op) :
      PDE_HelperDiscretization(global_op),
      static_matrices_initialized_(false) {
    Init_(plist);
  }

  PDE_Abstract(Teuchos::ParameterList& plist, Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      PDE_HelperDiscretization(mesh),
      static_matrices_initialized_(false) {
    global_op_ = Teuchos::null;
    Init_(plist);
  }
  ~PDE_Abstract() {};

  // main members 
  // -- required by the interface
  using PDE_HelperDiscretization::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p) override;
  // -- new interface for pre-computed data  
  void UpdateMatrices(double t);

  // -- setup
  void SetupTensor(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K) { K_ = K; }
  void SetupPoly(const Teuchos::RCP<const std::vector<WhetStone::Polynomial> >& K) { Kpoly_ = K; }
  void SetupPolyVector(const Teuchos::RCP<const std::vector<WhetStone::VectorPolynomial> >& K) { Kvec_ = K; }
  void Setup(const Teuchos::RCP<const std::vector<WhetStone::VectorSpaceTimePolynomial> >& K) {
    Kvec_st_ = K;
    if (!static_matrices_initialized_) CreateStaticMatrices_();
  }

  // optional calculation of flux from potential p
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override {};

 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;
  Teuchos::RCP<const std::vector<WhetStone::Polynomial> > Kpoly_;
  Teuchos::RCP<const std::vector<WhetStone::VectorPolynomial> > Kvec_;
  Teuchos::RCP<const std::vector<WhetStone::VectorSpaceTimePolynomial> > Kvec_st_;

  Schema global_schema_row_, global_schema_col_;
  Schema local_schema_col_, local_schema_row_;

 private:
  void Init_(Teuchos::ParameterList& plist);
  void CreateStaticMatrices_();

 private:
  std::string matrix_;
  bool grad_on_test_;

  Teuchos::RCP<WhetStone::BilinearForm> mfd_;

  bool static_matrices_initialized_;
  std::vector<std::vector<WhetStone::DenseMatrix> > static_matrices_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

