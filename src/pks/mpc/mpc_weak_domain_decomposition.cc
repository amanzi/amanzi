/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Multi process coupler for sequential coupling across subdomains.

#include "mpc_weak_domain_decomposition.hh"

namespace Amanzi {

// PK methods
// -- dt is the minimum of the sub pks
double
MPCWeakDomainDecomposition::get_dt()
{
  double sub_dt = 1.0e99;
  for (auto & pk : sub_pks_) {
    sub_dt = std::min<double>(sub_dt, pk->get_dt());
  }

  double dt;
  comm_->MinAll(&sub_dt, &dt, 1);
  return dt;
}

// Set timestep for sub PKs
void
MPCWeakDomainDecomposition::set_dt( double dt)
{
  for (auto& pk : sub_pks_) {
    pk->set_dt(dt);
  }
}

// -- advance each sub pk dt.
bool
MPCWeakDomainDecomposition::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool sub_fail = false;
  for (auto & pk : sub_pks_) {
    sub_fail = pk->AdvanceStep(t_old, t_new, reinit);
    if (sub_fail) break;
  }

  int sub_fail_i = sub_fail;
  int fail;
  comm_->MinAll(&sub_fail_i, &fail, 1);
  return (bool) fail;
}

// init sub PKs
void
MPCWeakDomainDecomposition::init_(const Teuchos::RCP<State>& S)
{
  // grab the list of subpks
  auto subpks = plist_->get<Teuchos::Array<std::string> >("PKs order");
  if (subpks.size() != 1) {
    Errors::Message msg;
    msg << "MPCWeakSubgrid: \"PKs order\" should consist of a single domain set of sub-pks.";
    Exceptions::amanzi_throw(msg);
  }
  std::string subgrid_name = subpks[0];
  subpks.pop_back();

  KeyTriple subgrid_triple;
  bool is_ds = Keys::splitDomainSet(subgrid_name, subgrid_triple);
  if (!is_ds) {
    Errors::Message msg;
    msg << "MPCWeakSubgrid: subpk \"" << subgrid_name << "\" should be a domain-set PK of the form SUBGRID_DOMAIN_NAME:*-NAME";
    Exceptions::amanzi_throw(msg);
  }

  // get the domain set and save the comm of the parent mesh for later
  const auto& ds = S->GetDomainSet(std::get<0>(subgrid_triple));
  comm_ = ds->parent->get_comm();

  // add one PK for each subdomain
  for (const auto& dname_id : *ds) {
    subpks.push_back(Keys::getKey(dname_id.first, std::get<2>(subgrid_triple)));
  }

  // -- create the lifted PKs
  PKFactory pk_factory;
  for (auto subpk : subpks) {
    // create the solution vector
    Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector());
    solution_->PushBack(pk_soln);

    // create the PK
    sub_pks_.push_back(pk_factory.CreatePK(subpk, pk_tree_, global_list_, S, pk_soln));
  }
}


} // namespace
