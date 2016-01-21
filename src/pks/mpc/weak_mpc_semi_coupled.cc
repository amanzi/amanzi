#include "Teuchos_XMLParameterListHelpers.hpp"

#include "pk_physical_bdf_base.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "strong_mpc.hh"

#include "weak_mpc_semi_coupled.hh"


namespace Amanzi {

// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------


/*

bool WeakMPCSemiCoupled::advance(double dt) {
  bool fail = false;
  MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();

  fail = (*pk)->advance(dt);


  const Epetra_MultiVector& surf_p = *S_next_->GetFieldData("surface-pressure")->ViewComponent("cell", false);
  Epetra_MultiVector& surfstar_p = *S_inter_->GetFieldData("surface_star-pressure",sub_pks_[1]->name())->ViewComponent("cell", false);
  

  for (unsigned c=0; c<surf_p.MyLength(); c++){
    if(surf_p[0][c] > 101325.0)
      surfstar_p[0][c] = surf_p[0][c];
    else
      surfstar_p[0][c]=101325.0;
  }

/--
  *S_inter_->GetFieldData("surface_star-pressure",sub_pks_[1]->name()) =
  *S_next_->GetFieldData("surface-pressure"); 
--/
  Teuchos::RCP<PKBDFBase> pk_domain =
    Teuchos::rcp_dynamic_cast<PKBDFBase>(sub_pks_[1]);
  ASSERT(pk_domain.get()); // make sure the pk_domain is not empty
  pk_domain->ChangedSolution(S_inter_.ptr());

if(fail) return fail;

  ++pk;

 fail = (*pk)->advance(dt);

  
  Epetra_MultiVector& surf_pr = *S_next_->GetFieldData("surface-pressure", S_inter_->GetField("surface-pressure")->owner())->ViewComponent("cell", false);
  const Epetra_MultiVector& surfstar_pr = *S_next_->GetFieldData("surface_star-pressure")->ViewComponent("cell", false);

  for (unsigned c=0; c<surf_pr.MyLength(); c++){
    if(surfstar_pr[0][c] > 101325.0)
      surf_pr[0][c] = surfstar_pr[0][c];
    else {}
    
  }
  
/--
  *S_next_->GetFieldData("surface-pressure",S_inter_->GetField("surface-pressure")->owner()) =
  *S_next_->GetFieldData("surface_star-pressure"); 
*--/
  CopySurfaceToSubsurface(*S_next_->GetFieldData("surface-pressure"),
			  S_next_->GetFieldData("pressure", S_inter_->GetField("pressure")->owner()).ptr());
  
  
  Teuchos::RCP<PKBDFBase> pk_surf =
    Teuchos::rcp_dynamic_cast<PKBDFBase>(sub_pks_[0]);
  ASSERT(pk_surf.get());
  pk_surf->ChangedSolution();
  
  return fail;
};


*/


//-------------------------------------------------------------------------------------
/*
// water coupling only -- merge later
bool WeakMPCSemiCoupled::advance(double dt) {
  bool fail = false;
  MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();

  // advance surface_star-pressure from t_n to t_(n+1)
  fail = (*pk)->advance(dt);
  
  Epetra_MultiVector& surf_pr = *S_inter_->GetFieldData("surface-pressure", S_inter_->GetField("surface-pressure")->owner())->ViewComponent("cell", false);
  const Epetra_MultiVector& surfstar_pr = *S_next_->GetFieldData("surface_star-pressure")->ViewComponent("cell", false);

  for (unsigned c=0; c<surf_pr.MyLength(); c++){
    if(surfstar_pr[0][c] > 101325.0)
      surf_pr[0][c] = surfstar_pr[0][c];
    else {}
  }
  
  CopySurfaceToSubsurface(*S_inter_->GetFieldData("surface-pressure"),
			  S_inter_->GetFieldData("pressure", S_inter_->GetField("pressure")->owner()).ptr());
  
  // NOTE: later do it in the setup --aj
  Teuchos::RCP<PKBDFBase> pk_domain =
    Teuchos::rcp_dynamic_cast<PKBDFBase>(sub_pks_[1]);
  ASSERT(pk_domain.get()); // make sure the pk_domain is not empty
  pk_domain->ChangedSolution(S_inter_.ptr());
  
  if(fail) return fail;  
  // advance surface-pressure from t_n to t_(n+1)
  ++pk;
  fail = (*pk)->advance(dt);
  if (fail) return fail;
  
  
  const Epetra_MultiVector& surf_p = *S_next_->GetFieldData("surface-pressure")->ViewComponent("cell", false);
  Epetra_MultiVector& surfstar_p = *S_next_->GetFieldData("surface_star-pressure",sub_pks_[0]->name())->ViewComponent("cell", false);
  
  for (unsigned c=0; c<surf_p.MyLength(); c++){
    if(surf_p[0][c] > 101325.0)
      surfstar_p[0][c] = surf_p[0][c];
    else
      surfstar_p[0][c]=101325.0;

    if(surf_p[0][c] > 101380.0)
      std::cout<<"PRESSURE-VALUE2: "<<surf_p[0][c]<<"\n";
  }
  
  // Mark surface_star-pressure evaluator as changed.
  // NOTE: later do it in the setup --aj
  Teuchos::RCP<PKBDFBase> pk_surf =
    Teuchos::rcp_dynamic_cast<PKBDFBase>(sub_pks_[0]);
  ASSERT(pk_surf.get());
  pk_surf->ChangedSolution();
  
  return fail;
};
 
*/





//-------------------------------------------------------------------------------------
// Semi coupled thermal hydrology
bool WeakMPCSemiCoupled::advance(double dt) {
  bool fail = false;
  MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();

  // advance surface_star-pressure from t_n to t_(n+1)
  fail = (*pk)->advance(dt);
  
  Epetra_MultiVector& surf_pr = *S_inter_->GetFieldData("surface-pressure", S_inter_->GetField("surface-pressure")->owner())->ViewComponent("cell", false);
  const Epetra_MultiVector& surfstar_pr = *S_next_->GetFieldData("surface_star-pressure")->ViewComponent("cell", false);

  for (unsigned c=0; c<surf_pr.MyLength(); c++){
    if(surfstar_pr[0][c] > 101325.0)
      surf_pr[0][c] = surfstar_pr[0][c];
    else {}
  }
  
  
  
  // copy surface_star-pressure at t_(n+1) to surface-pressure at t_n
  *S_inter_->GetFieldData("surface-temperature",S_inter_->GetField("surface-temperature")->owner()) =
  *S_next_->GetFieldData("surface_star-temperature"); 
  
  
  CopySurfaceToSubsurface(*S_inter_->GetFieldData("surface-pressure"),
			  S_inter_->GetFieldData("pressure", S_inter_->GetField("pressure")->owner()).ptr());

  CopySurfaceToSubsurface(*S_inter_->GetFieldData("surface-temperature"),
			 S_inter_->GetFieldData("temperature", S_inter_->GetField("temperature")->owner()).ptr());
  
  // NOTE: later do it in the setup --aj
  Teuchos::RCP<PKBDFBase> pk_domain =
    Teuchos::rcp_dynamic_cast<PKBDFBase>(sub_pks_[1]);
  ASSERT(pk_domain.get()); // make sure the pk_domain is not empty
  pk_domain->ChangedSolution(S_inter_.ptr());
  
  if(fail) return fail;  
  // advance surface-pressure from t_n to t_(n+1)
  ++pk;
  fail += (*pk)->advance(dt);
  if (fail) return fail;
  
  
  const Epetra_MultiVector& surf_p = *S_next_->GetFieldData("surface-pressure")->ViewComponent("cell", false);
  Epetra_MultiVector& surfstar_p = *S_next_->GetFieldData("surface_star-pressure",
				   S_inter_->GetField("surface_star-pressure")->owner())->ViewComponent("cell", false);
  
  for (unsigned c=0; c<surf_p.MyLength(); c++){
    if(surf_p[0][c] > 101325.0)
      surfstar_p[0][c] = surf_p[0][c];
    else
      surfstar_p[0][c]=101325.0;
  }
  
  *S_next_->GetFieldData("surface_star-temperature",S_inter_->GetField("surface_star-temperature")->owner()) = 
    *S_next_->GetFieldData("surface-temperature");


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


