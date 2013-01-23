/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for flow and energy.  This couples using a
block-diagonal coupler.
------------------------------------------------------------------------- */

#ifndef PKS_MPC_DIAGONAL_FLOW_ENERGY_HH_
#define PKS_MPC_DIAGONAL_FLOW_ENERGY_HH_

#include "richards.hh"
#include "two_phase.hh"
#include "strong_mpc.hh"

namespace Amanzi {

class MPCDiagonalFlowEnergy : public StrongMPC {

 public:
  MPCDiagonalFlowEnergy(Teuchos::ParameterList& plist,
                        const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      StrongMPC(plist, soln) {}

  // Virtual destructor
  virtual ~MPCDiagonalFlowEnergy() {}

  // initialize the preconditioner
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // advance
  virtual bool advance(double dt);

  // evaluate the residual
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

 protected: // methods
  void precon_diagonal(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  void precon_upper_triangular(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  void precon_lower_triangular(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  void precon_alternating(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  void precon_accumulation(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

 protected: // data

  Teuchos::RCP<Flow::Richards> richards_pk_;
  Teuchos::RCP<Energy::TwoPhase> two_phase_pk_;
  Teuchos::RCP<Epetra_MultiVector> D_pT_;
  Teuchos::RCP<Epetra_MultiVector> D_Tp_;

  enum PreconMethod {PRECON_ACCUMULATION,
                     PRECON_DIAGONAL,
                     PRECON_UPPER_TRIANGULAR,
                     PRECON_LOWER_TRIANGULAR,
                     PRECON_ALTERNATING};

  PreconMethod method_;
  double damping_;
  int n_iter_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCDiagonalFlowEnergy> reg_;

};

} // namespace


#endif
