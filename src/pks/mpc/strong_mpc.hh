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

#include "bdf_fn_base.hh"
#include "mpc.hh"
#include "pk_bdf_base.hh"

namespace Amanzi {

// note this looks odd, but StrongMPC is both a MPC within a hierarchy of BDF
// PKs, but it also IS a BDF PK itself, in that it implements the BDF
// interface.
class StrongMPC : public MPC<PKBDFBase>, public PKBDFBase {

public:
  StrongMPC(Teuchos::ParameterList& plist,
            const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      MPC<PKBDFBase>(plist,soln),
      PKBDFBase(plist,soln) {}

  virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // StrongMPC is a BDFFnBase
  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk fun().
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // -- enorm for the coupled system
  virtual double enorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du);

  // StrongMPC's preconditioner is, by default, just the block-diagonal
  // operator formed by placing the sub PK's preconditioners on the diagonal.
  // -- Apply preconditioner to u and returns the result in Pu.
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // -- Update the preconditioner.
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

  // -- Experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in
  //    state.
  virtual void changed_solution();

  // -- Admissibility of the solution.
  virtual bool is_admissible(Teuchos::RCP<const TreeVector> u);

private:
  // factory registration
  static RegisteredPKFactory<StrongMPC> reg_;

};
} // close namespace Amanzi

#endif
