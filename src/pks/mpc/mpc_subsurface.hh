/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for coupling energy and water in the subsurface,
with freezing.

------------------------------------------------------------------------- */

#ifndef MPC_SUBSURFACE_HH_
#define MPC_SUBSURFACE_HH_

#include "mpc_coupled_cells.hh"

namespace Amanzi {

class MPCDelegateEWCSubsurface;

class MPCSubsurface : public MPCCoupledCells {

 public:

  MPCSubsurface(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                Teuchos::ParameterList& FElist,
                const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, FElist, soln),
      MPCCoupledCells(plist, FElist, soln) {}

  // -- Initialize owned (dependent) variables.
  virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void initialize(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);

  // update the predictor to be physically consistent
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> up0,
          Teuchos::RCP<TreeVector> up);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // preconditioner application
  virtual void ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du);

 protected:

  enum PreconditionerType {
    PRECON_NONE = 0,
    PRECON_BLOCK_DIAGONAL = 1,
    PRECON_PICARD = 2,
    PRECON_EWC = 3
  };

  // preconditioner methods
  PreconditionerType precon_type_;

  // EWC delegate
  Teuchos::RCP<MPCDelegateEWCSubsurface> ewc_;

  bool dumped_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCSubsurface> reg_;

};


} // namespace

#endif


