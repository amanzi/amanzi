/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Riemann flux-based advection operator for a scalar field.
  The family will include DG and CG methods for linear and nonlinear
  fluxes. We shall start with the existing node2-face1 schema.
*/

#ifndef AMANZI_OPERATOR_PDE_ADVECTION_RIEMANN_HH_
#define AMANZI_OPERATOR_PDE_ADVECTION_RIEMANN_HH_

#include <string>

#include "VectorPolynomial.hh"

#include "PDE_Advection.hh"

namespace Amanzi {
namespace Operators {

class PDE_AdvectionRiemann : public PDE_Advection {
 public:
  PDE_AdvectionRiemann(Teuchos::ParameterList& plist,
                       Teuchos::RCP<Operator> global_op) :
      PDE_Advection(plist, global_op)
  {
    InitAdvection_(plist);
  }

  PDE_AdvectionRiemann(Teuchos::ParameterList& plist,
                       Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      PDE_Advection(plist, mesh)
  {
    InitAdvection_(plist);
  }

  // main members 
  // -- setup
  virtual void Setup(const CompositeVector& u) {};

  void Setup(const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> >& Kc,
             const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& Kf) {
    Kc_ = Kc;
    Kf_ = Kf;
  }

  // -- generate linearized operator: standard interface
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p) {};

  // -- generate linearized operator: new interface
  void UpdateMatrices(const Teuchos::Ptr<const std::vector<WhetStone::Polynomial> >& u);

  // -- determine advected flux of potential u
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& h,
                          const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::RCP<BCs>& bc,
                          const Teuchos::Ptr<CompositeVector>& flux);

 private:
  void InitAdvection_(Teuchos::ParameterList& plist);

 private:
  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > Kc_;
  Teuchos::RCP<std::vector<WhetStone::Polynomial> > Kf_;

  std::string method_, matrix_, flux_;
  int method_order_;
  bool jump_on_test_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

