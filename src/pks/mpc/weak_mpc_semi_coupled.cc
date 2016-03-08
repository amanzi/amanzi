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



/*

//-------------------------------------------------------------------------------------
// Semi coupled thermal hydrology -- FULL 3D RUN
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
 
*/
void
WeakMPCSemiCoupled::setup(const Teuchos::Ptr<State>& S) {
  S->AliasMesh("surface", "surface_star");
  MPC<PK>::setup(S);


  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string> >("PKs order");
  numPKs_ = names.size();  
  
  coupling_key_ = plist_->get<std::string>("coupling key"," ");
  ASSERT(!(coupling_key_ == " "));
};




//-------------------------------------------------------------------------------------
// Semi coupled thermal hydrology
bool WeakMPCSemiCoupled::advance(double dt) {
  bool fail = false;
 
  if (coupling_key_ == "surface subsurface system: columns"){
    fail = CoupledSurfSubsurfColumns(dt);
  }
  else if(coupling_key_ == "surface subsurface system: 3D"){
    fail = CoupledSurfSubsurf3D(dt);
  }
  
  return fail;
  
};
 
bool
WeakMPCSemiCoupled::CoupledSurfSubsurfColumns(double dt){
  bool fail = false;
  MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
    
    // advance surface_star-pressure from t_n to t_(n+1)
    fail = (*pk)->advance(dt);
    
    //copying surface_star (2D) data (pressures/temperatures) to column surface (1D-cells)[all the surf column cells get updates]
    const Epetra_MultiVector& surfstar_pres = *S_next_->GetFieldData("surface_star-pressure")->ViewComponent("cell", false);
    const Epetra_MultiVector& surfstar_temp = *S_next_->GetFieldData("surface_star-temperature")->ViewComponent("cell", false);
    
    unsigned int size_t = surfstar_pres.MyLength();
    
    //copying pressure
    for (unsigned c=0; c<size_t; c++){
      if(surfstar_pres[0][c] > 101325.0){
	std::stringstream name;
	name << "column_" << c <<"_surface";
	Epetra_MultiVector& surf_pres = *S_inter_->GetFieldData(getKey(name.str(),"pressure"), 
								S_inter_->GetField(getKey(name.str(),"pressure"))->owner())->ViewComponent("cell", false);
	surf_pres[0][0] = surfstar_pres[0][c];
      }
      else {}
    }
    
    
    //copying temperatures
    for (unsigned c=0; c<size_t; c++){
      std::stringstream name, name_ss;
      name << "column_" << c <<"_surface";
      name_ss << "column_" << c;
      Epetra_MultiVector& surf_temp = *S_inter_->GetFieldData(getKey(name.str(),"temperature"), 
							      S_inter_->GetField(getKey(name.str(),"temperature"))->owner())->ViewComponent("cell", false);
      surf_temp[0][0] = surfstar_temp[0][c];
      
      
      CopySurfaceToSubsurface(*S_inter_->GetFieldData(getKey(name.str(),"pressure")),
			      S_inter_->GetFieldData(getKey(name_ss.str(),"pressure"), 
						     S_inter_->GetField(getKey(name_ss.str(),"pressure"))->owner()).ptr());
      
      CopySurfaceToSubsurface(*S_inter_->GetFieldData(getKey(name.str(),"temperature")),
			      S_inter_->GetFieldData(getKey(name_ss.str(),"temperature"), 
						     S_inter_->GetField(getKey(name_ss.str(),"temperature"))->owner()).ptr());
    } 
    // NOTE: later do it in the setup --aj
    
    
    for(int i=1; i<numPKs_; i++){
      Teuchos::RCP<PKBDFBase> pk_domain =
	Teuchos::rcp_dynamic_cast<PKBDFBase>(sub_pks_[i]);
      ASSERT(pk_domain.get()); // make sure the pk_domain is not empty
      pk_domain->ChangedSolution(S_inter_.ptr());
    }
    if(fail) return fail;  
    
    // advance surface-pressure from t_n to t_(n+1)
    ++pk;
    while (pk != sub_pks_.end()){
      fail += (*pk)->advance(dt);
      if (fail) return fail;
      ++pk;
    }
    
    Epetra_MultiVector& surfstar_p = *S_next_->GetFieldData("surface_star-pressure",
							    S_inter_->GetField("surface_star-pressure")->owner())->ViewComponent("cell", false);
    Epetra_MultiVector& surfstar_t = *S_next_->GetFieldData("surface_star-temperature",
							    S_inter_->GetField("surface_star-temperature")->owner())->ViewComponent("cell", false);
    
    for (unsigned c=0; c<size_t; c++){
      std::stringstream name;
      name << "column_" << c <<"_surface";
      const Epetra_MultiVector& surf_p = *S_next_->GetFieldData(getKey(name.str(),"pressure"))->ViewComponent("cell", false);
      if(surf_p[0][c] > 101325.0)
	surfstar_p[0][c] = surf_p[0][0];
      else
	surfstar_p[0][c]=101325.0;
    }
    
    for (unsigned c=0; c<size_t; c++){
      std::stringstream name;
      name << "column_" << c <<"_surface";
      const Epetra_MultiVector& surf_t = *S_next_->GetFieldData(getKey(name.str(),"temperature"))->ViewComponent("cell", false);
      surfstar_t[0][c] = surf_t[0][0];
    }
    
    
    // Mark surface_star-pressure evaluator as changed.
    // NOTE: later do it in the setup --aj
    Teuchos::RCP<PKBDFBase> pk_surf =
      Teuchos::rcp_dynamic_cast<PKBDFBase>(sub_pks_[0]);
    ASSERT(pk_surf.get());
    pk_surf->ChangedSolution();
    return fail;
}


bool WeakMPCSemiCoupled::CoupledSurfSubsurf3D(double dt) {
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







} // namespace Amanzi


