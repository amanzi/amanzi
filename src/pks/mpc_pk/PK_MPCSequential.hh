/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Konstantin Lipnikov

  Implementation for the derived PK_MPCSequential class. Provides only the 
  AdvanceStep() method missing from MPC.hh. In sequential coupling, we 
  iteratively loop over the sub-PKs, calling their AdvanceStep() methods 
  until a strong convergence achieved and returning failure if any fail.

  See additional documentation in the base class src/pks/mpc_pk/PK_MPC.hh
*/

#ifndef AMANZI_PK_MPC_SEQUENTIAL_HH_
#define AMANZI_PK_MPC_SEQUENTIAL_HH_

#include "PK_MPC.hh"
#include "PK.hh"

namespace Amanzi {

class PK_MPCSequential : public PK_MPC<PK> {
 public:
  PK_MPCSequential(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt){};

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

  // access
  int num_itrs() { return num_itrs_; }
  double error_norm() { return error_norm_; }

 protected:
  int max_itrs_, num_itrs_;
  double error_norm_, tol_;
};

} // namespace Amanzi

#endif
