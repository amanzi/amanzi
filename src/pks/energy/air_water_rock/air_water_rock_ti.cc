/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "air_water_rock.hh"

namespace Amanzi {
namespace Energy {

// AirWaterRock is a BDFFnBase
// computes the non-linear functional f = f(t,u,udot)
void AirWaterRock::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                 Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) {
  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);

  // pointer-copy temperature into states and update any auxilary data
  solution_to_state(u_old, S_inter_);
  solution_to_state(u_new, S_next_);
  UpdateSecondaryVariables_();

  // get access to the solution
  Teuchos::RCP<CompositeVector> res = f->data();

  // conduction term

  // source term?

  // accumulation term
  AddAccumulation_(res);

  // advection term
  AddAdvection_(res);
};

// applies preconditioner to u and returns the result in Pu
void AirWaterRock::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  *Pu = *u;
};

// computes a norm on u-du and returns the result
double AirWaterRock::enorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du) {
  double enorm_val = 0.0;
  Teuchos::RCP<const Epetra_MultiVector> temp_vec = u->data()->ViewComponent("cell", false);
  Teuchos::RCP<const Epetra_MultiVector> dtemp_vec = du->data()->ViewComponent("cell", false);

  for (unsigned int lcv=0; lcv != temp_vec->MyLength(); ++lcv) {
    double tmp = abs((*(*dtemp_vec)(0))[lcv])/(atol_ + rtol_*abs((*(*temp_vec)(0))[lcv]));
    enorm_val = std::max<double>(enorm_val, tmp);
  }

#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  return enorm_val;
};

// updates the preconditioner
void AirWaterRock::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {}

} // namespace Energy
} // namespace Amanzi
