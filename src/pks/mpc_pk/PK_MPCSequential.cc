/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Konstantin Lipnikov
*/

/*
  MPC PK

  Implementation for the derived PK_MPCSequential class. Provides only the
  AdvanceStep() method missing from MPC.hh. In sequential coupling, we
  iteratively loop over the sub-PKs, calling their AdvanceStep() methods
  until a strong convergence achieved and returning failure if any fail.

  See additional documentation in the base class src/pks/mpc_pk/PK_MPC.hh
*/

#include "Key.hh"
#include "LScheme_Helpers.hh"

#include "PK_MPCSequential.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

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
  tol_ = sublist->get<double>("error tolerance", 1e-6);

  // pass L-scheme options to sub-pks
  L_scheme_ = sublist->get<bool>("L-scheme stabilization", false);
  if (L_scheme_) {
    L_scheme_keys_ = SetupLSchemeKey();
  }

  vo_ = Teuchos::rcp(new VerboseObject(soln->Comm(), "MPC_Sequential", *global_list));
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
// Inform PKs about the failure of the iterative coupling.
// -----------------------------------------------------------------------------
void
PK_MPCSequential::set_dt(double dt)
{
  for (auto pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) (*pk)->set_dt(dt);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
PK_MPCSequential::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail(false);
  TreeVector solution_copy(*solution_);

  // initialize L-scheme
  InitializeLSchemeStep();
  ComputeLSchemeStability();

  num_itrs_ = 0;
  error_norm_ = 0.0;

  while (num_itrs_ < 2 || (error_norm_ > tol_ && num_itrs_ < max_itrs_)) {
    auto du = Teuchos::rcp(new TreeVector(*solution_));

    for (int i = 0; i < sub_pks_.size(); ++i) {
      fail = sub_pks_[i]->AdvanceStep(t_old, t_new, reinit);
      if (fail) {
        *solution_ = solution_copy;
        for (int k = 0; k < i; ++k) {
          sub_pks_[k]->FailStep(t_old, t_new, Tags::DEFAULT);
        }
        return fail;
      }
    }
    CommitSequentialStep(du, solution_);

    du->Update(-1.0, *solution_, 1.0);
    error_norm_ = ErrorNorm(solution_, du);
    num_itrs_++;

    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "sequential itrs #" << num_itrs_ << " error=" << error_norm_ << "\n";
    }
  }

  if (error_norm_ > tol_) {
    set_dt((t_new - t_old) / 2);
    fail = true;
  }

  return fail;
}


// -----------------------------------------------------------------------------
// Relative l2 norm is the default metric.
// -----------------------------------------------------------------------------
double
PK_MPCSequential::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
{
  double err(0.0), err_i, unorm_i, r;
  for (int i = 0; i < u->size(); ++i) { 
    du->SubVector(i)->Norm2(&err_i);
    u->SubVector(i)->Norm2(&unorm_i);
    err = std::max(err, err_i / unorm_i);
  }

  // Redefine error norm using the true residual
  if (L_scheme_) {
    err = 0.0;
    auto& data = S_->GetW<LSchemeData>("l_scheme_data", "state");

    for (auto& item : data) {
      item.second.seq_error[0] = item.second.last_step_residual;
      r = item.second.update();
      err = std::max(err, item.second.last_step_residual);

      item.second.print(std::cout);
      item.second.shift();
    }
  }
  return err;
}


// -----------------------------------------------------------------------------
// Update previous accumulation field
// -----------------------------------------------------------------------------
void
PK_MPCSequential::CommitSequentialStep(Teuchos::RCP<const TreeVector> u_old,
                                       Teuchos::RCP<const TreeVector> u_new)
{
  if (L_scheme_) {
    for (int i = 0; i < sub_pks_.size(); ++i) {
      if (!L_scheme_keys_[i].empty()) {
        Key key_prev = L_scheme_keys_[i] + "_prev";
        const auto& u1_c = *u_new->SubVector(i)->Data()->ViewComponent("cell");
        *S_->GetW<CV_t>(key_prev, "state").ViewComponent("cell") = u1_c;
      }
    }
  }
}

} // namespace Amanzi
