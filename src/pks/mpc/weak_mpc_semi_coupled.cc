#include "Teuchos_XMLParameterListHelpers.hpp"
#include "EpetraExt_RowMatrixOut.h"

#include "mpc_surface_subsurface_helpers.hh"
#include "weak_mpc_semi_coupled.hh"



namespace Amanzi {

/*
// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double WeakMPCSemiCoupled::get_dt() {
  double dt = 1.0e99;
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    dt = std::min<double>(dt, (*pk)->get_dt());
  }
  return dt;
};

*/
// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool WeakMPCSemiCoupled::advance(double dt) {
 
  bool fail = false;
  MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
  //advance surface_star-pressure from t_n to t_(n+1)
 
  fail = (*pk)->advance(dt);

  //--*S_inter_->GetFieldData("surface-pressure","surface flow") = *S_next_->GetFieldData("surface_star-pressure",sub_pks_[0]->name());
  // copy surface_star-pressure at t_(n+1) to surface-pressure at t_n

  //CopySurfaceToSubsurface(*S_inter_->GetFieldData("surface-pressure", "surface flow"),S_inter_->GetFieldData("pressure", "subsurface flow").ptr()); 

  
  ++pk;

 fail += (*pk)->advance(dt);


 *S_next_->GetFieldData("surface_star-pressure",sub_pks_[0]->name()) =  *S_next_->GetFieldData("surface-pressure","surface flow");

//later do it in the setup or constructor
 Teuchos::RCP<PKPhysicalBDFBase> pk_surf = Teuchos::rcp_dynamic_cast<PKPhysicalBDFBase>(sub_pks_[0]);
 pk_surf->ChangedSolution();
   
  if (fail)
    return fail;

};

 


void
WeakMPCSemiCoupled::setup(const Teuchos::Ptr<State>& S) {
  S->AliasMesh("surface", "surface_star");
  MPC<PK>::setup(S);
};

/*
void
WeakMPCSemiCoupled::initialize(const Teuchos::Ptr<State>& S) {
  MPC<PK>::initialize(S);
}
*/

} // namespace


