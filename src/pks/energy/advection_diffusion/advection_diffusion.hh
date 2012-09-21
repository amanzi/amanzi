/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Start of a process kernel for the energy equation to be used in thermal
permafrost.  This starts with the simplification that T > T_freezing, limiting
us to the air-water system.
------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_DIFFUSION_HH_
#define PKS_ENERGY_DIFFUSION_HH_

#include "pk_factory.hh"
#include "bdf_time_integrator.hh"

#include "advection.hh"
#include "matrix_mfd.hh"

#include "pk_physical_bdf_base.hh"

namespace Amanzi {
namespace Energy {

class AdvectionDiffusion : public PKPhysicalBDFBase {

public:

  AdvectionDiffusion(Teuchos::ParameterList& plist,
                     const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist,solution),
      PKPhysicalBDFBase(plist, solution) {}

  // AdvectionDiffusion is a PK
  // -- Setup data.
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {}

  // -- Calculate any diagnostics prior to doing vis
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}


  // AdvectionDiffusion is a BDFFnBase
  // computes the non-linear functional f = f(t,u,udot)
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

private:
  // helper methods for calling the above methods
  void AddAccumulation_(const Teuchos::RCP<CompositeVector> f);
  void AddAdvection_(const Teuchos::RCP<State> S,
                     const Teuchos::RCP<CompositeVector> f, bool negate);
  void ApplyDiffusion_(const Teuchos::RCP<State> S, const Teuchos::RCP<CompositeVector> f);

  // methods for applying/using the discretization/operators
  void UpdateBoundaryConditions_();
  void ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature);

  // misc setup information
  Teuchos::ParameterList energy_plist_;
  double dt_;

  // boundary conditions
  Teuchos::RCP<Functions::BoundaryFunction> bc_temperature_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;
  std::vector<Operators::Matrix_bc> bc_markers_;
  std::vector<double> bc_values_;

  // operators
  Teuchos::RCP<Operators::Advection> advection_;
  Teuchos::RCP<Operators::MatrixMFD> matrix_;
  Teuchos::RCP<Operators::MatrixMFD> preconditioner_;

  // time integration
  Teuchos::RCP<BDFTimeIntegrator> time_stepper_;
  double atol_;
  double rtol_;

  // factory registration
  static RegisteredPKFactory<AdvectionDiffusion> reg_;

};

} // namespace Energy
} // namespace Amanzi

#endif
