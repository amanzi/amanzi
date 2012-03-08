/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "air_water_rock.hh"

namespace Amanzi {
namespace Energy {

// Energy_AirWaterRock is a BDFFnBase
// computes the non-linear functional f = f(t,u,udot)
virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                 Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) {
  S_inter_->set_time(t_old);
  S_->set_time(t_new);

  // pointer-copy temperature into states and update any auxilary data
  solution_to_state(u_old, S_inter_);
  solution_to_state(u_new, S_next_);
  UpdateSecondaryVariables_();

  // get access to the solution
  Teuchos::RCP<CompositeVector> res = f->get_data();

  // conduction term

  // source term?

  // accumulation term
  AddAccumulation_(f);

  // advection term
  AddAdvection_(f);
};

// applies preconditioner to u and returns the result in Pu
virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

// computes a norm on u-du and returns the result
virtual double enorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);

// updates the preconditioner
virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

} // namespace Energy
} // namespace Amanzi
