/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A generalized surface balance model using Crank-Nicholson time integration
for a system of ODEs, i.e.:

d theta(u)
---------  = Qin - Qout
  dt 


 ------------------------------------------------------------------------- */

#ifndef PK_SURFACE_BALANCE_BASE_HH_
#define PK_SURFACE_BALANCE_BASE_HH_

#include "pk_factory_ats.hh"
#include "pk_physical_bdf_base.hh"

namespace Amanzi {
namespace SurfaceBalance {

class SurfaceBalanceBase : public PKPhysicalBDFBase {

 public:
  SurfaceBalanceBase(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                         Teuchos::ParameterList& FElist,
                         const Teuchos::RCP<TreeVector>& solution);

  // main methods
  // -- Setup data.
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  virtual void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);
  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  
 protected:

  Key layer_;
  bool conserved_quantity_;
  Key source_key_;

  double theta_;

  Teuchos::RCP<CompositeVector> jac_;
  
 private:
  // factory registration
  static RegisteredPKFactory_ATS<SurfaceBalanceBase> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

