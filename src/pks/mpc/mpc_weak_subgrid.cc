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


#include "mpc_weak_subgrid.hh"


namespace Amanzi {


MPCWeakSubgrid::MPCWeakSubgrid(Teuchos::ParameterList& FElist,
        const Teuchos::RCP<Teuchos::ParameterList>& plist,
        const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution),
      MPC<PK>(FElist, plist, S, solution) {
  MPCWeakSubgrid::init_(S);
};


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double MPCWeakSubgrid::get_dt() {
  double dt = 1.0e99;
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    dt = std::min<double>(dt, (*pk)->get_dt());
  }
  return dt;
};

// -----------------------------------------------------------------------------
// Set timestep for sub PKs 
// -----------------------------------------------------------------------------
void MPCWeakSubgrid::set_dt( double dt) {
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_dt(dt);
  }

};

// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool MPCWeakSubgrid::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    fail = (*pk)->AdvanceStep(t_old, t_new, reinit);
    if (fail) {
      return fail;
    }
  }
  return fail;
};


void
MPCWeakSubgrid::init_(const Teuchos::RCP<State>& S)
{
  // grab the list of subpks
  auto subpks = plist_->get<Teuchos::Array<std::string> >("PKs order");
  AMANZI_ASSERT(subpks.size() == 1);
  std::string subgrid_name = subpks[0];
  subpks.pop_back();

  KeyTriple subgrid_triple;
  bool is_ds = Keys::splitDomainSet(subgrid_name, subgrid_triple);
  if (!is_ds) {
    Errors::Message msg;
    msg << "MPCWeakSubgrid: subpk \"" << subgrid_name << "\" should be a domain-set PK of the form SUBGRID_DOMAIN_NAME_*-NAME";
    Exceptions::amanzi_throw(msg);
  }

  // add for the various columns based on GIDs of the surface system
  std::string parent_domain = plist_->get<std::string>("parent domain");
  AmanziMesh::Entity_kind kind =
      AmanziMesh::entity_kind(plist_->get<std::string>("entity kind"));
  std::string regionname = plist_->get<std::string>("subgrid region name");
  
  Teuchos::RCP<const AmanziMesh::Mesh> parent_mesh = S->GetMesh(parent_domain);
  int nsubgrid = parent_mesh->num_entities(kind, AmanziMesh::Parallel_type::OWNED);
  auto& map = parent_mesh->map(kind,false);
  for (int i=0; i!=nsubgrid; ++i) {
    int gid = map.GID(i);
    std::stringstream domain_name_stream;
    domain_name_stream << std::get<0>(subgrid_triple) << "_" << gid;
    subpks.push_back(Keys::getKey(domain_name_stream.str(), std::get<2>(subgrid_triple)));
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
