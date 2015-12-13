#include "Teuchos_XMLParameterListHelpers.hpp"

#include "pk_physical_bdf_base.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "strong_mpc.hh"

#include "weak_mpc_semi_coupled.hh"


namespace Amanzi {

// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool WeakMPCSemiCoupled::advance(double dt) {
  bool fail = false;
  MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();

  // advance surface_star-pressure from t_n to t_(n+1)
  fail = (*pk)->advance(dt);
/*
  // copy surface_star-pressure at t_(n+1) to surface-pressure at t_n
  *S_inter_->GetFieldData("surface-pressure",S_inter_->GetField("surface-pressure")->owner()) =
    *S_next_->GetFieldData("surface_star-pressure");
  CopySurfaceToSubsurface(*S_inter_->GetFieldData("surface-pressure"),
			  S_inter_->GetFieldData("pressure", S_inter_->GetField("pressure")->owner()).ptr());
  // NOTE: later do it in the setup --aj
  Teuchos::RCP<PKBDFBase> pk_domain =
    Teuchos::rcp_dynamic_cast<PKBDFBase>(sub_pks_[1]);
  ASSERT(pk_domain.get()); // make sure the pk_domain is not empty
  pk_domain->ChangedSolution(S_inter_.ptr());
*/
if(fail) return fail;  
  // advance surface-pressure from t_n to t_(n+1)
  ++pk;
  fail += (*pk)->advance(dt);
  if (fail) return fail;
  
  // copy surface-pressure at t_(n+1) to surface_star-pressure at t_n+1
  *S_next_->GetFieldData("surface_star-pressure",sub_pks_[0]->name()) =
    *S_next_->GetFieldData("surface-pressure");
  
  // Mark surface_star-pressure evaluator as changed.
  // NOTE: later do it in the setup --aj
  Teuchos::RCP<PKBDFBase> pk_surf =
    Teuchos::rcp_dynamic_cast<PKBDFBase>(sub_pks_[0]);
  ASSERT(pk_surf.get());
  pk_surf->ChangedSolution();
  
  return fail;
};
 


void
WeakMPCSemiCoupled::setup(const Teuchos::Ptr<State>& S) {
  S->AliasMesh("surface", "surface_star");
  MPC<PK>::setup(S);
};


} // namespace Amanzi


