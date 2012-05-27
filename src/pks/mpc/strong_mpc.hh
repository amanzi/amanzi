/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived StrongMPC class.  Is both a PK and a Model
Evalulator, providing needed methods for BDF time integration of the coupled
system.

Completely automated and generic to any sub PKs, this uses a block diagonal
preconditioner.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#ifndef PKS_MPC_WEAKMPC_HH_
#define PKS_MPC_WEAKMPC_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "state.hh"
#include "tree_vector.hh"
#include "bdf_fn_base.hh"
#include "bdf_time_integrator.hh"
#include "mpc.hh"

namespace Amanzi {

class StrongMPC : public MPC {

public:
  StrongMPC(Teuchos::ParameterList& mpc_plist, const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& soln);

  virtual bool advance(double dt);
  virtual void initialize(const Teuchos::RCP<State>& S);
  virtual double get_dt() { return dt_; }

  // StrongMPC is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // computes a norm on u-du and returns the result
  virtual double enorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

protected:
  // mathematical operators
  Teuchos::RCP<Amanzi::BDFTimeIntegrator> time_stepper_;
  double dt_;
  double atol_;
  double rtol_;
  double time_step_reduction_factor_;

private:
  // factory registration
  static RegisteredPKFactory<StrongMPC> reg_;

};
} // close namespace Amanzi

#endif
