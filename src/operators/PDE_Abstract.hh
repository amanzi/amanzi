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

#include "Teuchos_RCPStdSharedPtrConversions.hpp"

// Amanzi
#include "BilinearForm.hh"
#include "CoefficientModel.hh"
#include "InterfaceWhetStone.hh"
#include "Polynomial.hh"
#include "Tensor.hh"
#include "VectorObjects.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "OperatorDefs.hh"
#include "PDE_HelperDiscretization.hh"

namespace Amanzi {
namespace Operators {

class PDE_Abstract : public PDE_HelperDiscretization {
 public:
  PDE_Abstract(Teuchos::ParameterList& plist, Teuchos::RCP<Operator> global_op) :
      PDE_HelperDiscretization(global_op),
      coef_type_(CoefType::CONSTANT),
      static_matrices_initialized_(false) {
    Init_(plist);
  }

  PDE_Abstract(Teuchos::ParameterList& plist, Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      PDE_HelperDiscretization(mesh),
      coef_type_(CoefType::CONSTANT),
      static_matrices_initialized_(false) {
    global_op_ = Teuchos::null;
    Init_(plist);
  }
  ~PDE_Abstract() {};

  // main members 
  // -- required by the interface
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
          const Teuchos::Ptr<const CompositeVector>& p=Teuchos::null);
  void UpdateMatrices() {
    UpdateMatrices(Teuchos::null, Teuchos::null);
  }

  // -- new interface for pre-computed data  
  void UpdateMatrices(double t);

  // -- setup can be used to change coefficient type before any call 
  //    of UpdateMatrices. Note that pointers to previous coefficient
  //    values are not deleted
  template<typename T>
  void Setup(const Teuchos::RCP<std::vector<T> >& K, bool reset);

 protected:
  // available models for operator coefficient
  Teuchos::RCP<std::vector<WhetStone::Polynomial> > Kpoly_;
  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > Kvec_poly_;
  Teuchos::RCP<std::vector<WhetStone::VectorSpaceTimePolynomial> > Kvec_stpoly_;

  Schema global_schema_row_, global_schema_col_;
  Schema local_schema_col_, local_schema_row_;

 private:
  void Init_(Teuchos::ParameterList& plist);
  void CreateStaticMatrices_();

 private:
  std::string matrix_;
  bool grad_on_test_;

  Teuchos::RCP<WhetStone::BilinearForm> mfd_;
  Teuchos::RCP<InterfaceWhetStone> interface_;

  CoefType coef_type_;
  bool static_matrices_initialized_;
  std::vector<std::vector<WhetStone::DenseMatrix> > static_matrices_;
};


/* ******************************************************************
* Specialization of Setup
****************************************************************** */
template<>
inline
void PDE_Abstract::Setup<WhetStone::Tensor>(
    const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K, bool reset) {
  coef_type_ = CoefType::CONSTANT;

  auto Kc = Teuchos::get_shared_ptr(K);
  const auto coef = std::make_shared<CoefficientModel<WhetStone::Tensor> >(Kc);
  interface_ = Teuchos::rcp(new InterfaceWhetStoneMFD<
      WhetStone::BilinearForm, CoefficientModel<WhetStone::Tensor> >(mfd_, coef));
}
 
template<>
inline
void PDE_Abstract::Setup<WhetStone::Polynomial>(
    const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& K, bool reset) {
  Kpoly_ = K;
  coef_type_ = CoefType::POLYNOMIAL;

  auto Kc = Teuchos::get_shared_ptr(K);
  const auto coef = std::make_shared<CoefficientModel<WhetStone::Polynomial> >(Kc);
  interface_ = Teuchos::rcp(new InterfaceWhetStoneMFD<
      WhetStone::BilinearForm, CoefficientModel<WhetStone::Polynomial> >(mfd_, coef));
}

template<>
inline
void PDE_Abstract::Setup<WhetStone::VectorPolynomial>(
    const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> >& K, bool reset) {
  Kvec_poly_ = K;
  coef_type_ = CoefType::VECTOR_POLYNOMIAL;
}

template<>
inline
void PDE_Abstract::Setup<WhetStone::VectorSpaceTimePolynomial>(
    const Teuchos::RCP<std::vector<WhetStone::VectorSpaceTimePolynomial> >& K, bool reset) {
  Kvec_stpoly_ = K;
  coef_type_ = CoefType::VECTOR_SPACETIME_POLYNOMIAL;
  if (!static_matrices_initialized_ || reset) CreateStaticMatrices_();
}

}  // namespace Operators
}  // namespace Amanzi

#endif

