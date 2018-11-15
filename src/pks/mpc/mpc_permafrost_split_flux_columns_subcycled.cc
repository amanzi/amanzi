/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

------------------------------------------------------------------------- */

#include "primary_variable_field_evaluator.hh"
#include "mpc_surface_subsurface_helpers.hh"

#include "mpc_permafrost_split_flux_columns_subcycled.hh"

#include "PK_Physical.hh"

namespace Amanzi {

MPCPermafrostSplitFluxColumnsSubcycled::MPCPermafrostSplitFluxColumnsSubcycled(Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution),
      MPC<PK>(FElist, plist, S, solution)
{
  // collect keys and names
  std::string domain = plist_->get<std::string>("domain name");
  std::string domain_star = plist_->get<std::string>("star domain name");

  // generate a list of sub-pks, based upon the star system + a list of columns
  // -- get the original PKs order -- the column must be the last one
  auto subpks = plist_->get<Teuchos::Array<std::string> >("PKs order");
  std::string colname = subpks[subpks.size()-1];
  subpks.pop_back();

  // -- generate the column triple domain set
  KeyTriple col_triple;
  bool is_ds = Keys::splitDomainSet(colname, col_triple);
  if (!is_ds) {
    Errors::Message msg;
    msg << "WeakMPCSemiCoupled subpk: \"" << colname << "\" should be a domain-set PK of the form column_*-NAME";
    Exceptions::amanzi_throw(msg);
  }

  // -- add for the various columns based on GIDs of the surface system
  auto surf_mesh = S->GetMesh(domain_star);
  int ncols = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int i=0; i!=ncols; ++i) {
    int gid = surf_mesh->cell_map(false).GID(i);
    std::stringstream domain_name_stream;
    domain_name_stream << std::get<0>(col_triple) << "_" << gid;
    subpks.push_back(Keys::getKey(domain_name_stream.str(), std::get<2>(col_triple)));

    std::stringstream surf_domain_name_stream;
    surf_domain_name_stream << domain << "_" << gid;
    col_domains_.push_back(surf_domain_name_stream.str());
  }

  // set up keys
  p_primary_variable_suffix_ = plist_->get<std::string>("pressure primary variable suffix", "pressure");
  T_primary_variable_suffix_ = plist_->get<std::string>("temperature primary variable suffix", "temperature");

  p_primary_variable_star_ = Keys::readKey(*plist_, domain_star, "pressure primary variable star", Keys::getVarName(p_primary_variable_suffix_));
  T_primary_variable_star_ = Keys::readKey(*plist_, domain_star, "temperature primary variable star", Keys::getVarName(T_primary_variable_suffix_));

  p_lateral_flow_source_suffix_ = plist_->get<std::string>("mass lateral flow source suffix", "mass_lateral_flow_source");
  T_lateral_flow_source_suffix_ = plist_->get<std::string>("energy lateral flow source suffix", "energy_lateral_flow_source");

  p_conserved_variable_star_ = Keys::readKey(*plist_, domain_star, "mass conserved quantity star", "water_content");
  T_conserved_variable_star_ = Keys::readKey(*plist_, domain_star, "energy conserved quantity star", "energy");

  cv_key_ = Keys::readKey(*plist_, domain_star, "cell volume", "cell_volume");
  
  // set up for a primary variable field evaluator for the flux
  Key primary_pkey = Keys::getKey(domain + std::string("_*"), p_lateral_flow_source_suffix_);
  std::cout << "setting eval for " << primary_pkey << std::endl;
  auto& p_sublist = S->FEList().sublist(primary_pkey);
  p_sublist.set("field evaluator type", "primary variable");
  Key primary_Tkey = Keys::getKey(domain + std::string("_*"), T_lateral_flow_source_suffix_);
  std::cout << "setting eval for " << primary_Tkey << std::endl;
  auto& T_sublist = S->FEList().sublist(primary_Tkey);
  T_sublist.set("field evaluator type", "primary variable");

  // init sub-pks
  plist_->set("PKs order", subpks);
  init_(S);
};


void MPCPermafrostSplitFluxColumnsSubcycled::Initialize(const Teuchos::Ptr<State>& S)
{
  // initialize the columns
  for (int i=1; i!=sub_pks_.size(); ++i) sub_pks_[i]->Initialize(S);

  // copy the columns to the star system, which initializes the star system consistently without an IC
  CopyPrimaryToStar(S, S);
  S->GetField(p_primary_variable_star_, S->GetField(p_primary_variable_star_)->owner())->set_initialized();
  S->GetField(T_primary_variable_star_, S->GetField(T_primary_variable_star_)->owner())->set_initialized();

  // initialize the star system (IC already set, but other things might happen)
  sub_pks_[0]->Initialize(S);
}


void MPCPermafrostSplitFluxColumnsSubcycled::Setup(const Teuchos::Ptr<State>& S)
{
  MPC<PK>::Setup(S);

  // require the coupling sources
  for (const auto& col_domain : col_domains_) {
    Key pkey = Keys::getKey(col_domain, p_lateral_flow_source_suffix_);
    S->RequireField(pkey, pkey)
        ->SetMesh(S->GetMesh(col_domain))
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(pkey);

    Key Tkey = Keys::getKey(col_domain, T_lateral_flow_source_suffix_);
    S->RequireField(Tkey, Tkey)
        ->SetMesh(S->GetMesh(col_domain))
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(Tkey);
  }
}

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double MPCPermafrostSplitFluxColumnsSubcycled::get_dt()
{
  return sub_pks_[0]->get_dt();
};

// -----------------------------------------------------------------------------
// Set timestep for sub PKs 
// -----------------------------------------------------------------------------
void MPCPermafrostSplitFluxColumnsSubcycled::set_dt( double dt)
{
  AMANZI_ASSERT(false);
  sub_pks_[0]->set_dt(dt);
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_dt(dt);
  }
};

// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool MPCPermafrostSplitFluxColumnsSubcycled::AdvanceStep(double t_old, double t_new, bool reinit)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // Advance the star system 
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Beginning timestepping on star system." << std::endl;

  bool fail = false;
  fail = sub_pks_[0]->AdvanceStep(t_old, t_new, reinit);
  fail |= !sub_pks_[0]->ValidStep();
  if (fail) return fail;
  sub_pks_[0]->CommitStep(t_old, t_new, S_next_);

  // Copy star's new value into primary's old value
  CopyStarToPrimary(t_new - t_old);

  // Now advance the primary
  for (int i=1; i!=sub_pks_.size(); ++i) {
    double t_inner = t_old;
    bool done = false;
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "Beginning timestepping on column " << i << std::endl;

    S_inter_->set_time(t_old);
    while (!done) {
      double dt_inner = std::min(sub_pks_[i]->get_dt(), t_new - t_inner);
      *S_next_->GetScalarData("dt", "coordinator") = dt_inner;
      S_next_->set_time(t_inner + dt_inner);
      bool fail_inner = sub_pks_[i]->AdvanceStep(t_inner, t_inner+dt_inner, false);
      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "  step failed? " << fail_inner << std::endl;
      fail_inner |= !sub_pks_[i]->ValidStep();
      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "  step failed or was not valid? " << fail_inner << std::endl;

      if (fail_inner) {
        dt_inner = sub_pks_[i]->get_dt();
        *S_next_ = *S_inter_;

        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "  failed, new timestep is " << dt_inner << std::endl;

        
      } else {
        sub_pks_[i]->CommitStep(t_inner, t_inner + dt_inner, S_next_);
        t_inner += dt_inner;
        *S_inter_ = *S_next_;
        if (t_inner >= t_new - 1.e-10) {
          done = true;
        }
        dt_inner = sub_pks_[i]->get_dt();
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "  success, new timestep is " << dt_inner << std::endl;
      }
    }
  }
  S_inter_->set_time(t_old);

  // Copy the primary into the star to advance
  CopyPrimaryToStar(S_next_.ptr(), S_next_.ptr());

  return false;
};

bool MPCPermafrostSplitFluxColumnsSubcycled::ValidStep() 
{
  return true;
}
  


void MPCPermafrostSplitFluxColumnsSubcycled::CommitStep(double t_old, double t_new,
        const Teuchos::RCP<State>& S)
{}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFluxColumnsSubcycled::CopyPrimaryToStar(const Teuchos::Ptr<const State>& S,
                                    const Teuchos::Ptr<State>& S_star)
{
  // copy p primary variables into star primary variable
  auto& p_star = *S_star->GetFieldData(p_primary_variable_star_, S_star->GetField(p_primary_variable_star_)->owner())
                  ->ViewComponent("cell",false);
  for (int c=0; c!=p_star.MyLength(); ++c) {
    Key pkey = Keys::getKey(col_domains_[c], p_primary_variable_suffix_);
    const auto& p = *S->GetFieldData(pkey)->ViewComponent("cell",false);
    AMANZI_ASSERT(p.MyLength() == 1);
    if (p[0][0] <= 101325.0) {
      p_star[0][c] = 101325.;
    } else {
      p_star[0][c] = p[0][0];
    }
  }

  auto peval = S_star->GetFieldEvaluator(p_primary_variable_star_);
  auto peval_pvfe = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(peval);
  peval_pvfe->SetFieldAsChanged(S_star.ptr());

  // copy T primary variable
  auto& T_star = *S_star->GetFieldData(T_primary_variable_star_, S_star->GetField(T_primary_variable_star_)->owner())
                  ->ViewComponent("cell",false);
  for (int c=0; c!=p_star.MyLength(); ++c) {
    Key Tkey = Keys::getKey(col_domains_[c], T_primary_variable_suffix_);
    const auto& T = *S->GetFieldData(Tkey)->ViewComponent("cell",false);
    AMANZI_ASSERT(T.MyLength() == 1);
    T_star[0][c] = T[0][0];
  }

  auto Teval = S_star->GetFieldEvaluator(T_primary_variable_star_);
  auto Teval_pvfe = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(Teval);
  Teval_pvfe->SetFieldAsChanged(S_star.ptr());
}

// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFluxColumnsSubcycled::CopyStarToPrimary(double dt)
{
  // make sure we have the evaluator at the new state timestep
  if (p_eval_pvfes_.size() == 0) {
    for (const auto& col_domain : col_domains_) {
      Key pkey = Keys::getKey(col_domain, p_lateral_flow_source_suffix_);
      Teuchos::RCP<FieldEvaluator> fe = S_next_->GetFieldEvaluator(pkey);
      p_eval_pvfes_.push_back(Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe));
      AMANZI_ASSERT(p_eval_pvfes_.back() != Teuchos::null);
    }
  }
  
  if (T_eval_pvfes_.size() == 0) {
    for (const auto& col_domain : col_domains_) {
      Key pkey = Keys::getKey(col_domain, T_lateral_flow_source_suffix_);
      Teuchos::RCP<FieldEvaluator> fe = S_next_->GetFieldEvaluator(pkey);
      T_eval_pvfes_.push_back(Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe));
      AMANZI_ASSERT(T_eval_pvfes_.back() != Teuchos::null);
    }
  }

  // these updates should do nothing, but you never know
  S_inter_->GetFieldEvaluator(p_conserved_variable_star_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetFieldEvaluator(p_conserved_variable_star_)->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator(T_conserved_variable_star_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetFieldEvaluator(T_conserved_variable_star_)->HasFieldChanged(S_next_.ptr(), name_);

  // grab the data, difference
  Epetra_MultiVector q_div(*S_next_->GetFieldData(p_conserved_variable_star_)->ViewComponent("cell",false));
  q_div.Update(1.0/dt,
               *S_next_->GetFieldData(p_conserved_variable_star_)->ViewComponent("cell",false),
               -1.0/dt,
               *S_inter_->GetFieldData(p_conserved_variable_star_)->ViewComponent("cell",false),
               0.);

  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(1.0, *S_next_->GetFieldData(cv_key_)->ViewComponent("cell",false), q_div, 0.);

  // copy into columns
  for (int c=0; c!=q_div.MyLength(); ++c) {
    Key pkey = Keys::getKey(col_domains_[c], p_lateral_flow_source_suffix_);
    (*S_next_->GetFieldData(pkey, pkey)->ViewComponent("cell",false))[0][0] = q_div[0][c];
    p_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());
  }
  
  // grab the data, difference
  Epetra_MultiVector qE_div(*S_next_->GetFieldData(T_conserved_variable_star_)->ViewComponent("cell",false));
  qE_div.Update(1.0/dt,
               *S_next_->GetFieldData(T_conserved_variable_star_)->ViewComponent("cell",false),
               -1.0/dt,
               *S_inter_->GetFieldData(T_conserved_variable_star_)->ViewComponent("cell",false),
               0.);

  // scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(1.0, *S_next_->GetFieldData(cv_key_)->ViewComponent("cell",false), qE_div, 0.);

  // copy into columns
  for (int c=0; c!=qE_div.MyLength(); ++c) {
    Key Tkey = Keys::getKey(col_domains_[c], T_lateral_flow_source_suffix_);
    (*S_next_->GetFieldData(Tkey, Tkey)->ViewComponent("cell",false))[0][0] = qE_div[0][c];
    T_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());
  }
}

} // namespace
