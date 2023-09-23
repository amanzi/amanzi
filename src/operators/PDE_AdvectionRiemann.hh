/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Riemann flux-based advection operator for a scalar field.
  The family will include DG and CG methods for linear and nonlinear
  fluxes. We shall start with the existing node2-face1 schema.
*/

#ifndef AMANZI_OPERATOR_PDE_ADVECTION_RIEMANN_HH_
#define AMANZI_OPERATOR_PDE_ADVECTION_RIEMANN_HH_

#include <string>

#include "DG_Modal.hh"
#include "SpaceTimePolynomial.hh"
#include "VectorObjects.hh"

#include "PDE_Advection.hh"

namespace Amanzi {
namespace Operators {

struct SurfaceFluxData {
  WhetStone::DenseMatrix Uface, Dface;
  double uflux, dflux;
};


class PDE_AdvectionRiemann : public PDE_Advection {
 public:
  PDE_AdvectionRiemann(Teuchos::ParameterList& plist, Teuchos::RCP<Operator> global_op)
    : PDE_Advection(plist, global_op), static_matrices_initialized_(false)
  {
    InitAdvection_(plist);
  }

  PDE_AdvectionRiemann(Teuchos::ParameterList& plist, Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : PDE_Advection(plist, mesh), static_matrices_initialized_(false)
  {
    InitAdvection_(plist);
  }

  // main members
  // -- setup for various flux algorithms
  virtual void Setup(const CompositeVector& u) final{};

  void Setup(const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial>>& Kc,
             const Teuchos::RCP<std::vector<WhetStone::Polynomial>>& Kf)
  {
    Kc_ = Kc;
    Kf_ = Kf;
  }
  void Setup(const Teuchos::Ptr<const std::vector<WhetStone::SpaceTimePolynomial>>& uc, bool reset)
  {
    uc_ = uc;
    if (!static_matrices_initialized_ || reset) CreateStaticMatrices_();
  }

  // -- generate linearized operator: standard interface
  using PDE_Advection::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p) final{};

  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              double (*)(double)) final{};

  // -- generate linearized operator: new interface
  void UpdateMatrices(const std::vector<WhetStone::Polynomial>& u);
  void UpdateMatrices(double t);

  // -- determine advected flux of potential u
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& h,
                          const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::RCP<BCs>& bc,
                          const Teuchos::Ptr<CompositeVector>& flux) final;

  // boundary conditions
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) final;

  // access
  const WhetStone::DG_Modal& dg() const { return *dg_; }

 private:
  void InitAdvection_(Teuchos::ParameterList& plist);
  void CreateStaticMatrices_();

 private:
  std::string method_, matrix_, flux_;
  bool jump_on_test_;

  Teuchos::RCP<WhetStone::DG_Modal> dg_;

  // support of Rusanov flux
  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial>> Kc_;
  Teuchos::RCP<std::vector<WhetStone::Polynomial>> Kf_;

  // support of space-time polynomial velocity
  bool static_matrices_initialized_;
  Teuchos::Ptr<const std::vector<WhetStone::SpaceTimePolynomial>> uc_;
  std::vector<std::vector<SurfaceFluxData>> static_matrices_;
};

} // namespace Operators
} // namespace Amanzi

#endif
