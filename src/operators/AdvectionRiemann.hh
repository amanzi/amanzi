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
  // -- data
  virtual void UpdateMatrices(const CompositeVector& u);
  virtual void UpdateMatrices(const CompositeVector& u, const CompositeVector& dhdT) {};
  // -- results -- determine advected flux of u
  void UpdateFlux(const CompositeVector& h , const CompositeVector& u,
                  const Teuchos::RCP<BCs>& bc, CompositeVector& flux);
  
  // boundary conditions
  void ApplyBCs(bool primary, bool eliminate);

 private:
  void InitAdvection_(Teuchos::ParameterList& plist);
  void UpdateMatricesCell_(const CompositeVector& u);
  void UpdateMatricesFace_(const CompositeVector& u);

 private:
  std::string flux_, riemann_;

  // AdvectionFunction adv_velocity_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

