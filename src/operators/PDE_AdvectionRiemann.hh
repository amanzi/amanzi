/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_PDE_ADVECTION_RIEMANN_HH_
#define AMANZI_OPERATOR_PDE_ADVECTION_RIEMANN_HH_

#include <string>

#include "DG_Modal.hh"
#include "VectorPolynomial.hh"

#include "PDE_Advection.hh"

namespace Amanzi {
namespace Operators {

class PDE_AdvectionRiemann : public PDE_Advection {
 public:
  PDE_AdvectionRiemann(Teuchos::ParameterList& plist,
                       Teuchos::RCP<Operator> global_op)
    : PDE_Advection(plist, global_op)
  {
    InitAdvection_(plist);
  }

  PDE_AdvectionRiemann(Teuchos::ParameterList& plist,
                       Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : PDE_Advection(plist, mesh)
  {
    InitAdvection_(plist);
  }

  // main members
  // -- setup
  virtual void Setup(const CompositeVector& u) override{};

  void Setup(const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial>>& Kc,
             const Teuchos::RCP<std::vector<WhetStone::Polynomial>>& Kf)
  {
    Kc_ = Kc;
    Kf_ = Kf;
  }

  // -- generate linearized operator: standard interface
  virtual void
  UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                 const Teuchos::Ptr<const CompositeVector>& p) override{};

  // -- generate linearized operator: new interface
  void UpdateMatrices(
    const Teuchos::Ptr<const std::vector<WhetStone::Polynomial>>& u);

  // -- determine advected flux of potential u
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& h,
                          const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::RCP<BCs>& bc,
                          const Teuchos::Ptr<CompositeVector>& flux) override;

  // boundary conditions
  virtual void
  ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;

  // access
  const WhetStone::DG_Modal& dg() const { return *dg_; }

 private:
  void InitAdvection_(Teuchos::ParameterList& plist);

 private:
  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial>> Kc_;
  Teuchos::RCP<std::vector<WhetStone::Polynomial>> Kf_;

  std::string method_, matrix_, flux_;
  bool jump_on_test_;

  Teuchos::RCP<WhetStone::DG_Modal> dg_;
};

} // namespace Operators
} // namespace Amanzi

#endif
