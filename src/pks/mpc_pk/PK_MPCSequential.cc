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

#include "PK_MPCSequential.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
PK_MPCSequential::PK_MPCSequential(Teuchos::ParameterList& pk_tree,
                                   const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                   const Teuchos::RCP<State>& S,
                                   const Teuchos::RCP<TreeVector>& soln)
  : PK_MPC<PK>(pk_tree, global_list, S, soln)
{
  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  auto pk_list = Teuchos::sublist(global_list, "PKs", true);
  auto sublist = Teuchos::sublist(pk_list, pk_name, true);

  max_itrs_ = sublist->get<int>("maximum number of iterations", 100);
  tol_ = sublist->get<double>("error tolerance", 1e-6);
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double
PK_MPCSequential::get_dt()
{
  double dt = 1.0e99;
  for (PK_MPC<PK>::SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    dt = std::min<double>(dt, (*pk)->get_dt());
  }
  return dt;
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
PK_MPCSequential::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;

  // create copy of the solution
  TreeVector solution_copy(*solution_);

  // iterations require to reset the primary field, e.g. for correct
  // calculation of the accumulation term
  num_itrs_ = 0;
  error_norm_ = 1.0e+20;
  double sol_norm(1.0);

  while (error_norm_ > tol_ && num_itrs_ < max_itrs_) {
    TreeVector solution_tmp(*solution_);

    int i = 0;
    for (PK_MPC<PK>::SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
      *solution_->SubVector(i)->Data() = *solution_copy.SubVector(i)->Data();

      fail = (*pk)->AdvanceStep(t_old, t_new, reinit);
      if (fail) return fail;

      ++i;
    }

    // calculate error
    if (num_itrs_ > 0) {
      solution_tmp.Update(-1.0, *solution_, 1.0);
      solution_tmp.Norm2(&error_norm_);
      solution_->Norm2(&sol_norm);
      error_norm_ /= sol_norm;
    }
    num_itrs_++;
  }

  return fail;
}

} // namespace Amanzi
