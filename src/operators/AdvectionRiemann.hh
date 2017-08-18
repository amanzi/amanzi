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

#ifndef AMANZI_OPERATOR_ADVECTION_RIEMANN_HH_
#define AMANZI_OPERATOR_ADVECTION_RIEMANN_HH_

#include <string>

#include "Advection.hh"

namespace Amanzi {
namespace Operators {

class AdvectionRiemann : public Advection {
 public:
  AdvectionRiemann(Teuchos::ParameterList& plist,
                   Teuchos::RCP<Operator> global_op) :
      Advection(plist, global_op)
  {
    InitAdvection_(plist);
  }

  AdvectionRiemann(Teuchos::ParameterList& plist,
                   Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      Advection(plist, mesh)
  {
    InitAdvection_(plist);
  }

  // required members 
  // -- setup
  virtual void Setup(const CompositeVector& u) {};

  // -- generate linearized operator: standard interface
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p) {};

  // -- generate linearized operator: new interface interface
  void UpdateMatrices(const Teuchos::Ptr<const std::vector<WhetStone::Polynomial> >& u);

  // -- determine advected flux of potential u
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& h,
                          const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::RCP<BCs>& bc,
                          Teuchos::Ptr<CompositeVector>& flux);

 private:
  void InitAdvection_(Teuchos::ParameterList& plist);

 private:
  std::string flux_, riemann_;
  int method_order_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

