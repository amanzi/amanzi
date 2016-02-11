/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Interface for the derived MPCWeak class.  Provides only the Advance() method
  missing from MPC.hh.  In weak coupling, we simply loop over the sub-PKs,
  calling their Advance() methods and returning failure if any fail.

  Simplest form of sequential coupling.

  See additional documentation in the base class src/pks/mpc_pk/MPC_PK.hh
*/

#ifndef AMANZI_WEAK_MPC_HH_
#define AMANZI_WEAK_MPC_HH_

#include "MPC_PK.hh"
#include "PK.hh"

namespace Amanzi {

class MPCWeak : public MPC_PK<PK> {
 public:
  MPCWeak(Teuchos::ParameterList& pk_tree,
          const Teuchos::RCP<Teuchos::ParameterList>& global_list,
          const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& soln) :
      PK_Default(pk_tree, global_list, S, soln),
      MPC_PK<PK>(pk_tree, global_list, S, soln) {};

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt){};

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false);

 private:
  // factory registration
  static RegisteredPKFactory<MPCWeak> reg_;
};

}  // namespace Amanzi

#endif
