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

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"

#include "state.hh"
#include "strong_mpc.hh"

namespace Amanzi {

class MPCDiagonalFlowEnergy : public StrongMPC {

 public:
  MPCDiagonalFlowEnergy(Teuchos::ParameterList& mpc_plist, const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& soln);

  // initialize the preconditioner
  virtual void initialize(const Teuchos::RCP<State>& S);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> up);

 protected:
  void precon_diagonal(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  void precon_upper_triangular(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  void precon_lower_triangular(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  void precon_alternating(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  void precon_accumulation(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  Teuchos::RCP<Epetra_MultiVector> D_pT_;
  Teuchos::RCP<Epetra_MultiVector> D_Tp_;

  enum PreconMethod {PRECON_ACCUMULATION,
                     PRECON_DIAGONAL,
                     PRECON_UPPER_TRIANGULAR,
                     PRECON_LOWER_TRIANGULAR,
                     PRECON_ALTERNATING};

  PreconMethod method_;
  double damping_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCDiagonalFlowEnergy> reg_;

};

} // namespace


#endif
