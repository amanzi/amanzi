/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived OperatorsplitMPC class.  Provides only the advance() method
missing from MPC.hh.  In operatorsplit coupling, we simply loop over the sub-PKs,
calling their advance() methods and returning failure if any fail.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#ifndef PKS_MPC_OPERATORSPLITMPC_HH_
#define PKS_MPC_OPERATORSPLITMPC_HH_

#include "PK.hh"
#include "mpc.hh"

namespace Amanzi {

class OperatorSplitMPC : public MPC<PK> {

 public:

  OperatorSplitMPC(Teuchos::ParameterList& FElist,
          const Teuchos::RCP<Teuchos::ParameterList>& plist,
          const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~OperatorSplitMPC() = default;

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();

  // -- initialize in reverse order
  virtual void Initialize(const Teuchos::Ptr<State>& S) {
    sub_pks_[1]->Initialize(S);
    sub_pks_[0]->Initialize(S);
  }
  
  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

  virtual void set_dt(double dt);

  virtual void CopyPrimaryToStar(const Teuchos::RCP<const State>& S,
          const Teuchos::RCP<State>& S_star);
  virtual void CopyStarToPrimary(const Teuchos::RCP<const State>& S_star,
          const Teuchos::RCP<State>& S);

 protected:
  Key primary_variable_;
  Key primary_variable_star_;
  
 private:
  // factory registration
  static RegisteredPKFactory<OperatorSplitMPC> reg_;


};
} // close namespace Amanzi

#endif
