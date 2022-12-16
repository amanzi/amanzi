/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  This is the mpc_pk component of the Amanzi code.

  Interface for the derived PK_MPCWeak class.  Provides only the Advance() method
  missing from MPC.hh.  In weak coupling, we simply loop over the sub-PKs,
  calling their Advance() methods and returning failure if any fail.

  Simplest form of sequential coupling.

  See additional documentation in the base class src/pks/mpc_pk/PK_MPC.hh
*/

#ifndef AMANZI_PK_MPC_WEAK_HH_
#define AMANZI_PK_MPC_WEAK_HH_

#include "PK_MPC.hh"
#include "PK.hh"

namespace Amanzi {

class PK_MPCWeak : public PK_MPC<PK> {
 public:
  PK_MPCWeak(Teuchos::ParameterList& pk_tree,
             const Teuchos::RCP<Teuchos::ParameterList>& global_list,
             const Teuchos::RCP<State>& S,
             const Teuchos::RCP<TreeVector>& soln)
    : PK_MPC<PK>(pk_tree, global_list, S, soln){};

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt){};

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

 private:
  // factory registration
  static RegisteredPKFactory<PK_MPCWeak> reg_;
};

} // namespace Amanzi

#endif
