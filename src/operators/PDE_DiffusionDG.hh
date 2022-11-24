/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  DiffusionDG implements the Diffusion interface using
  discontinuous polynomials.
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_DG_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_DG_HH_

#include <memory>
#include <string>

// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CoefficientModel.hh"
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "InterfaceWhetStone.hh"
#include "MatrixObjects.hh"
#include "Tensor.hh"
#include "WhetStoneFunction.hh"

// Operators
#include "PDE_HelperDiscretization.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class PDE_DiffusionDG : public PDE_HelperDiscretization {
 public:
  PDE_DiffusionDG(Teuchos::ParameterList& plist, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_HelperDiscretization(mesh), Kf_(NULL)
  {
    global_op_ = Teuchos::null;
    pde_type_ = PDE_DIFFUSION_DG;
    Init_(plist);
  }

  // main virtual members
  // -- setup
  template <typename T>
  void
  Setup(const std::shared_ptr<std::vector<T>>& Kc, const std::shared_ptr<std::vector<double>>& Kf);

  // -- creation of an operator
  using PDE_HelperDiscretization::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p) override;

  // -- modify local matrices due to boundary conditions
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;

  // -- postprocessing: calculated flux u from potential p
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;

  // access
  const WhetStone::DG_Modal& dg() const { return *dg_; }

 private:
  virtual void Init_(Teuchos::ParameterList& plist);

 private:
  std::string matrix_;
  int method_order_, numi_order_;

  // penalty coefficeint
  std::shared_ptr<std::vector<double>> Kf_;

  // other operators
  Teuchos::RCP<Op> jump_up_op_, jump_pu_op_, penalty_op_;

  Teuchos::RCP<WhetStone::DG_Modal> dg_;
  Teuchos::RCP<InterfaceWhetStone> interface_;
};


/* ******************************************************************
* Speciation of a member function
****************************************************************** */
template <>
inline void
PDE_DiffusionDG::Setup(const std::shared_ptr<std::vector<WhetStone::Tensor>>& Kc,
                       const std::shared_ptr<std::vector<double>>& Kf)
{
  Kf_ = Kf;

  auto coef = std::make_shared<CoefficientModel<WhetStone::Tensor>>(Kc);
  interface_ = Teuchos::rcp(
    new InterfaceWhetStoneDG<WhetStone::DG_Modal, CoefficientModel<WhetStone::Tensor>>(dg_, coef));
}

template <>
inline void
PDE_DiffusionDG::Setup(const std::shared_ptr<std::vector<WhetStone::WhetStoneFunction*>>& Kc,
                       const std::shared_ptr<std::vector<double>>& Kf)
{
  Kf_ = Kf;

  auto coef = std::make_shared<CoefficientModel<WhetStone::WhetStoneFunction*>>(Kc);
  interface_ = Teuchos::rcp(
    new InterfaceWhetStoneDG<WhetStone::DG_Modal, CoefficientModel<WhetStone::WhetStoneFunction*>>(
      dg_, coef));
}

template <>
inline void
PDE_DiffusionDG::Setup(const std::shared_ptr<std::vector<WhetStone::MatrixPolynomial>>& Kc,
                       const std::shared_ptr<std::vector<double>>& Kf)
{
  Kf_ = Kf;

  auto coef = std::make_shared<CoefficientModel<WhetStone::MatrixPolynomial>>(Kc);
  interface_ = Teuchos::rcp(
    new InterfaceWhetStoneDG<WhetStone::DG_Modal, CoefficientModel<WhetStone::MatrixPolynomial>>(
      dg_, coef));
}

} // namespace Operators
} // namespace Amanzi


#endif
