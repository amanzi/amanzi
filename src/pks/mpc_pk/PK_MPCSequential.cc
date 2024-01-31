/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Konstantin Lipnikov
*/

/*
  This is the mpc_pk component of the Amanzi code.

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

  auto tmp1 = Teuchos::sublist(global_list, "PKs", true);
  auto tmp2 = Teuchos::sublist(tmp1, pk_name, true);
  auto sublist = Teuchos::sublist(tmp2, "time integrator", true);

  max_itrs_ = sublist->get<int>("maximum number of iterations", 100);
  tol_ = sublist->get<double>("error tolerance", 1e-5);
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
  bool fail(false);

  num_itrs_ = 0;
  error_norm_ = 1.0e+20;

  while (error_norm_ > tol_ && num_itrs_ < max_itrs_) {
    auto du = Teuchos::rcp(new TreeVector(*solution_));

    for (auto pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
      fail = (*pk)->AdvanceStep(t_old, t_new, reinit);
      if (fail) return fail;
    }

    // calculate error
    if (num_itrs_ > 0) {
      du->Update(-1.0, *solution_, 1.0);
      error_norm_ = ErrorNorm(solution_, du);
    }
    num_itrs_++;

    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "sequential iteration #" << num_itrs_ << " error=" << error_norm_ << "\n";
    }
  }

  return fail;
}


// -----------------------------------------------------------------------------
// Relative l2 norm is the default metric.
// -----------------------------------------------------------------------------
double
PK_MPCSequential::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
{
  double err, unorm;
  du->Norm2(&err);
  u->Norm2(&unorm);
  return err / unorm;
}

} // namespace Amanzi
