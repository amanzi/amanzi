/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! Mixin for subgrid model MPCs with dynamic number of PKs

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

  NOTE currently this is just a weak MPC, but should be generalized as a
  mixin in the new rewrite of PKs.

 */

#ifndef ATS_PK_MPC_WEAK_SUBGRID_HH_
#define ATS_PK_MPC_WEAK_SUBGRID_HH_

#include "mpc.hh"

namespace Amanzi {

class MPCWeakSubgrid : public MPC<PK> {
 public:

  MPCWeakSubgrid(Teuchos::ParameterList& FElist,
          const Teuchos::RCP<Teuchos::ParameterList>& plist,
          const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& solution);

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

  virtual void set_dt(double dt);

  
 protected:
  void init_(const Teuchos::RCP<State>& S);

 private:
  // factory registration
  static RegisteredPKFactory<MPCWeakSubgrid> reg_;

};

} // namespace Amanzi

#endif
