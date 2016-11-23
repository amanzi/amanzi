#include "Teuchos_XMLParameterListHelpers.hpp"

//#include "pk_physical_bdf_base.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "strong_mpc.hh"

#include "weak_mpc_semi_coupled.hh"



namespace Amanzi {

  unsigned WeakMPCSemiCoupled::flag_star = 0;
  unsigned WeakMPCSemiCoupled::flag_star_surf = 0;
// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their valid_step() method
// -----------------------------------------------------------------------------
/*bool WeakMPCSemiCoupled::valid_step() {
  bool valid_local = MPC<PK>::valid_step();
  int valid_int_local = valid_local ? 1 : 0;
  int valid_int = 0;
  S_->GetMesh("surface")->get_comm()->MinAll(&valid_int_local, &valid_int, 1);
  return valid_int == 0 ? false : true;
};
*/

  /*
// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double WeakMPCSemiCoupled::get_dt() {
  double dt = 1.0e99;
  
  std::vector<double> time_steps;
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    time_steps.push_back((*pk)->get_dt());
    dt = std::min<double>(dt,(*pk)->get_dt());
  }
  min_dt_ = dt;
  --  double dt_col = 0.0;
  for (int i=1; i<time_steps.size(); i++){
    dt_col = std::min(dt_col, time_steps[i]);
      }
      min_dt_ = *min_element(time_steps.begin(), time_steps.end());
    //dt = std::min(time_steps[0],col_dt);
   
  if (time_steps[0] < sync_time_)
    dt = time_steps[0];
  else if  (sync_time_ < min_dt_) 
    dt =  min_dt_;
  else
    dt = sync_time_;
    
  surf_dt_ = dt;
  double dt_local = dt;
  S_->GetMesh("surface")->get_comm()->MinAll(&dt_local, &dt, 1);
  return dt;
};
  */

double WeakMPCSemiCoupled::get_dt() {
  double dt = 1.0e99;
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    dt = std::min<double>(dt,(*pk)->get_dt());
  }
  
  double dt_local = dt;
  S_->GetMesh("surface")->get_comm()->MinAll(&dt_local, &dt, 1);
  
  return dt;

}

// -----------------------------------------------------------------------------
// Set timestep for sub PKs 
// -----------------------------------------------------------------------------
void WeakMPCSemiCoupled::set_dt(double dt) {
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_dt(dt);
  }

};

// -----------------------------------------------------------------------------
// Set up each PK
// -----------------------------------------------------------------------------



void
WeakMPCSemiCoupled::Setup(const Teuchos::Ptr<State>& S) {
  S->AliasMesh("surface", "surface_star");
  MPC<PK>::Setup(S);


  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string> >("PKs order");
  numPKs_ = names.size();  
  
  coupling_key_ = plist_->get<std::string>("coupling key"," ");
  subcycle_key_ =  plist_->get<bool>("subcycle",false);
  //  sync_time_ = plist_->get<double>("sync time"); //provide default value later!!
  ASSERT(!(coupling_key_ == " "));
};

void 
WeakMPCSemiCoupled::Initialize(const Teuchos::Ptr<State>& S){

  MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
  ++pk;
  for (pk; pk!=sub_pks_.end(); ++pk){
    (*pk)->Initialize(S);
  }
  
  MPC<PK>::Initialize(S);

}
//-------------------------------------------------------------------------------------
// Semi coupled thermal hydrology
  bool WeakMPCSemiCoupled::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;
  
  if (coupling_key_ == "surface subsurface system: columns"){
    fail = CoupledSurfSubsurfColumns(t_old, t_new, reinit);
  }
  else if(coupling_key_ == "surface subsurface system: 3D"){
    fail = CoupledSurfSubsurf3D(t_old, t_new, reinit);
  }
  
  return fail;
  
};
 
bool
WeakMPCSemiCoupled::CoupledSurfSubsurfColumns(double t_old, double t_new, bool reinit){
  bool fail = false;
  MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
  
  // advance surface_star-pressure from t_n to t_(n+1)

  //ensure the star solution is marked as changed when the subsurface columns fail
  if (flag_star || flag_star_surf){
    flag_star = 0;
    flag_star_surf=0;
    Teuchos::RCP<PK_PhysicalBDF_Default> pk_sfstar =
      Teuchos::rcp_dynamic_cast<PK_PhysicalBDF_Default>(sub_pks_[0]);
    ASSERT(pk_sfstar.get());
    pk_sfstar->ChangedSolution();
  }
  
  fail = (*pk)->AdvanceStep(t_old, t_new, reinit);
  int nfailed_surf = 0;
  if (fail)
    nfailed_surf++;
  
  int nfailed_local_sf = nfailed_surf;
  S_->GetMesh("surface")->get_comm()->SumAll(&nfailed_local_sf, &nfailed_surf, 1);
 
  if(nfailed_surf > 0){
    flag_star_surf=1;
    return true;
  }
  //copying surface_star (2D) data (pressures/temperatures) to column surface (1D-cells)[all the surf column cells get updates]
  const Epetra_MultiVector& surfstar_pres = *S_next_->GetFieldData("surface_star-pressure")->ViewComponent("cell", false);
  const Epetra_MultiVector& surfstar_temp = *S_next_->GetFieldData("surface_star-temperature")->ViewComponent("cell", false);
  
  unsigned int size_t = surfstar_pres.MyLength();

  ASSERT(size_t == numPKs_ -1); // check if the subsurface columns are equal to the surface cells
  
  
  //copying pressure
  for (unsigned c=0; c<size_t; c++){
    if(surfstar_pres[0][c] > 101325.00){
      std::stringstream name;
      int id = S_->GetMesh("surface")->cell_map(false).GID(c);
     
      name << "column_" << id <<"_surface";
      Epetra_MultiVector& surf_pres = *S_inter_->GetFieldData(getKey(name.str(),"pressure"), 
							      S_inter_->GetField(getKey(name.str(),"pressure"))->owner())->ViewComponent("cell", false);
      surf_pres[0][0] = surfstar_pres[0][c];
    }
    else {}
  }
  
  
  //copying temperatures
  for (unsigned c=0; c<size_t; c++){
    std::stringstream name, name_ss;
    int id = S_->GetMesh("surface")->cell_map(false).GID(c);
    name << "column_" << id <<"_surface";
    name_ss << "column_" << id;
   
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
    Teuchos::RCP<PK_PhysicalBDF_Default> pk_domain =
      Teuchos::rcp_dynamic_cast<PK_PhysicalBDF_Default>(sub_pks_[i]);
    ASSERT(pk_domain.get()); // make sure the pk_domain is not empty
    pk_domain->ChangedSolution(S_inter_.ptr());
  }

  
 
 
  // advance surface-pressure from t_n to t_(n+1)
  ++pk;
  
  int nfailed = 0;
  int count=0;
  double t0 = S_inter_->time();
  double t1 = S_next_->time();

    
  for (pk; pk!=sub_pks_.end(); ++pk){

    if(!subcycle_key_){  
      bool c_fail = (*pk)->AdvanceStep(t_old, t_new, reinit);
      if (c_fail) nfailed++;
    }
    else
      {
      std::stringstream name, name_ss;
      int id = S_->GetMesh("surface")->cell_map(false).GID(count);
      name << "column_" << id <<"_surface";
      name_ss << "column_" << id;
      
      double loc_dt =0;//revisit dt;      
          
      bool done = false;
      double t = 0;
      bool cyc_flag  = false;
      //      S_inter_->set_time(t0+t);
      //S_next_->set_time(t0 + t + loc_dt);
      //lets put dt=0 for a while---revisit
      double dt =0;
      while(!done){
	bool fail_pk = false; //----revisit= (*pk)->advance(loc_dt);
	//fail_pk |= !(*pk)->valid_step();
		
	if(!fail_pk){

	  //--revisit (*pk)->commit_state(loc_dt, S_next_);
	  t = t + loc_dt;
	  double loc_dt_old = loc_dt;

	  double loc_dt_new = (*pk)->get_dt();
	  if (dt - (t+loc_dt_new) < 1.)
	    loc_dt_new = dt -t;
	  else if ( t + loc_dt_new < dt && (t+2*loc_dt_new > dt ) )
	    loc_dt_new = 0.5*(dt - t);

	  loc_dt = std::min<double>(loc_dt_new, dt - t);
	  
	  if (loc_dt <= 1.e-10){
	    done = true;
	    if (cyc_flag)
	      S_inter_->set_time(t0);
	    }
	  else{

	  *S_inter_->GetFieldData(getKey(name_ss.str(),"pressure"), S_inter_->GetField(getKey(name_ss.str(),"pressure"))->owner()) = 
	    *S_next_->GetFieldData(getKey(name_ss.str(),"pressure"));
					     
	  
	  *S_inter_->GetFieldData(getKey(name_ss.str(),"temperature"), S_inter_->GetField(getKey(name_ss.str(),"temperature"))->owner()) = 
	    *S_next_->GetFieldData(getKey(name_ss.str(),"temperature")); 
	  
	  *S_inter_->GetFieldData(getKey(name.str(),"pressure"), S_inter_->GetField(getKey(name.str(),"pressure"))->owner()) = 
	    *S_next_->GetFieldData(getKey(name.str(),"pressure"));
					     
	  
	  *S_inter_->GetFieldData(getKey(name.str(),"temperature"), S_inter_->GetField(getKey(name.str(),"temperature"))->owner()) = 
	    *S_next_->GetFieldData(getKey(name.str(),"temperature"));
	  
	  *S_inter_->GetFieldData(getKey(name.str(),"snow_depth"), S_inter_->GetField(getKey(name.str(),"snow_depth"))->owner()) = 
	    *S_next_->GetFieldData(getKey(name.str(),"snow_depth"));

	   *S_inter_->GetFieldData(getKey(name.str(),"conducted_energy_source"), 
				   S_inter_->GetField(getKey(name.str(),"conducted_energy_source"))->owner()) = 
	     *S_next_->GetFieldData(getKey(name.str(),"conducted_energy_source"));
	  
	   *S_inter_->GetFieldData(getKey(name.str(),"mass_source"), 
				   S_inter_->GetField(getKey(name.str(),"mass_source"))->owner()) = 
	     *S_next_->GetFieldData(getKey(name.str(),"mass_source"));

	   *S_inter_->GetFieldData(getKey(name_ss.str(),"mass_source"), 
				   S_inter_->GetField(getKey(name_ss.str(),"mass_source"))->owner()) = 
	     *S_next_->GetFieldData(getKey(name_ss.str(),"mass_source"));

	   *S_inter_->GetFieldData(getKey(name.str(),"mass_source_temperature"), 
				   S_inter_->GetField(getKey(name.str(),"mass_source_temperature"))->owner()) = 
	    *S_next_->GetFieldData(getKey(name.str(),"mass_source_temperature"));

	    *S_inter_->GetFieldData(getKey(name.str(),"snow_temperature"), 
				   S_inter_->GetField(getKey(name.str(),"snow_temperature"))->owner()) = 
	    *S_next_->GetFieldData(getKey(name.str(),"snow_temperature"));
	  	  
	    *S_inter_->GetFieldData(getKey(name.str(),"snow_age"), 
				   S_inter_->GetField(getKey(name.str(),"snow_age"))->owner()) = 
	    *S_next_->GetFieldData(getKey(name.str(),"snow_age"));
	    
	    *S_inter_->GetFieldData(getKey(name.str(),"snow_density"), 
				   S_inter_->GetField(getKey(name.str(),"snow_density"))->owner()) = 
	    *S_next_->GetFieldData(getKey(name.str(),"snow_density"));

	    *S_inter_->GetFieldData(getKey(name.str(),"stored_SWE"), 
				   S_inter_->GetField(getKey(name.str(),"stored_SWE"))->owner()) = 
	    *S_next_->GetFieldData(getKey(name.str(),"stored_SWE"));

	   *S_inter_->GetFieldData(getKey(name_ss.str(),"surface_subsurface_flux"), 
				  S_inter_->GetField(getKey(name_ss.str(),"surface_subsurface_flux"))->owner()) = 
	     *S_next_->GetFieldData(getKey(name_ss.str(),"surface_subsurface_flux"));
	  
	   *S_inter_->GetFieldData(getKey(name_ss.str(),"saturation_liquid"), 
				   S_inter_->GetField(getKey(name_ss.str(),"saturation_liquid"))->owner()) = 
	     *S_next_->GetFieldData(getKey(name_ss.str(),"saturation_liquid"));
	   
	   *S_inter_->GetFieldData(getKey(name.str(),"unfrozen_fraction"), 
				   S_inter_->GetField(getKey(name.str(),"unfrozen_fraction"))->owner()) = 
	     *S_next_->GetFieldData(getKey(name.str(),"unfrozen_fraction"));
	   
	   *S_inter_->GetFieldData(getKey(name.str(),"porosity"), 
				   S_inter_->GetField(getKey(name.str(),"porosity"))->owner()) = 
	     *S_next_->GetFieldData(getKey(name.str(),"porosity"));

	   *S_inter_->GetFieldData(getKey(name.str(),"evaporative_flux"), 
				   S_inter_->GetField(getKey(name.str(),"evaporative_flux"))->owner()) = 
	     *S_next_->GetFieldData(getKey(name.str(),"evaporative_flux"));
	   
	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe1 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S_inter_->GetFieldEvaluator(getKey(name.str(),"mass_source_temperature")));
	   
	   pfe1->SetFieldAsChanged(S_inter_.ptr());

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe2 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S_inter_->GetFieldEvaluator(getKey(name.str(),"conducted_energy_source")));
	   
	   pfe2->SetFieldAsChanged(S_inter_.ptr());

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe3 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S_inter_->GetFieldEvaluator(getKey(name.str(),"mass_source")));
	   
	   pfe3->SetFieldAsChanged(S_inter_.ptr());

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe4 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S_inter_->GetFieldEvaluator(getKey(name_ss.str(),"mass_source")));
	   
	   pfe4->SetFieldAsChanged(S_inter_.ptr());
	  
	   Teuchos::RCP<PK_PhysicalBDF_Default> pk_domain1 =
	     Teuchos::rcp_dynamic_cast<PK_PhysicalBDF_Default>(sub_pks_[count+1]);
	   ASSERT(pk_domain1.get()); // make sure the pk_domain is not empty
	   pk_domain1->ChangedSolution(S_inter_.ptr());
	   
	   S_inter_->set_time(t0+t);
	   S_next_->set_time( t0 + t + loc_dt);
	   
	  }
	  
	}
	else{
	  cyc_flag = true;
	  *S_next_->GetFieldData(getKey(name_ss.str(),"pressure"), S_next_->GetField(getKey(name_ss.str(),"pressure"))->owner()) = 
	    *S_inter_->GetFieldData(getKey(name_ss.str(),"pressure"));
					     
	  
	  *S_next_->GetFieldData(getKey(name_ss.str(),"temperature"), S_next_->GetField(getKey(name_ss.str(),"temperature"))->owner()) = 
	    *S_inter_->GetFieldData(getKey(name_ss.str(),"temperature")); 
	  
	  *S_next_->GetFieldData(getKey(name.str(),"pressure"), S_next_->GetField(getKey(name.str(),"pressure"))->owner()) = 
	    *S_inter_->GetFieldData(getKey(name.str(),"pressure"));
					     
	  
	  *S_next_->GetFieldData(getKey(name.str(),"temperature"), S_next_->GetField(getKey(name.str(),"temperature"))->owner()) = 
	    *S_inter_->GetFieldData(getKey(name.str(),"temperature"));
	  
	  *S_next_->GetFieldData(getKey(name.str(),"snow_depth"), S_next_->GetField(getKey(name.str(),"snow_depth"))->owner()) = 
	    *S_inter_->GetFieldData(getKey(name.str(),"snow_depth"));

	   *S_next_->GetFieldData(getKey(name.str(),"conducted_energy_source"), 
				   S_next_->GetField(getKey(name.str(),"conducted_energy_source"))->owner()) = 
	     *S_inter_->GetFieldData(getKey(name.str(),"conducted_energy_source"));
	  
	   *S_next_->GetFieldData(getKey(name.str(),"mass_source"), 
				   S_next_->GetField(getKey(name.str(),"mass_source"))->owner()) = 
	     *S_inter_->GetFieldData(getKey(name.str(),"mass_source"));

	   *S_next_->GetFieldData(getKey(name_ss.str(),"mass_source"), 
				   S_next_->GetField(getKey(name_ss.str(),"mass_source"))->owner()) = 
	     *S_inter_->GetFieldData(getKey(name_ss.str(),"mass_source"));

	   *S_next_->GetFieldData(getKey(name.str(),"mass_source_temperature"), 
				   S_next_->GetField(getKey(name.str(),"mass_source_temperature"))->owner()) = 
	    *S_inter_->GetFieldData(getKey(name.str(),"mass_source_temperature"));

	    *S_next_->GetFieldData(getKey(name.str(),"snow_temperature"), 
				   S_next_->GetField(getKey(name.str(),"snow_temperature"))->owner()) = 
	    *S_inter_->GetFieldData(getKey(name.str(),"snow_temperature"));
	  	  
	    *S_next_->GetFieldData(getKey(name.str(),"snow_age"), 
				   S_next_->GetField(getKey(name.str(),"snow_age"))->owner()) = 
	    *S_inter_->GetFieldData(getKey(name.str(),"snow_age"));
	    
	    *S_next_->GetFieldData(getKey(name.str(),"snow_density"), 
				   S_next_->GetField(getKey(name.str(),"snow_density"))->owner()) = 
	    *S_inter_->GetFieldData(getKey(name.str(),"snow_density"));

	    *S_next_->GetFieldData(getKey(name.str(),"stored_SWE"), 
				   S_next_->GetField(getKey(name.str(),"stored_SWE"))->owner()) = 
	    *S_inter_->GetFieldData(getKey(name.str(),"stored_SWE"));

	   *S_next_->GetFieldData(getKey(name_ss.str(),"surface_subsurface_flux"), 
				  S_next_->GetField(getKey(name_ss.str(),"surface_subsurface_flux"))->owner()) = 
	     *S_inter_->GetFieldData(getKey(name_ss.str(),"surface_subsurface_flux"));
	  
	  

	   *S_next_->GetFieldData(getKey(name_ss.str(),"saturation_liquid"), 
				  S_next_->GetField(getKey(name_ss.str(),"saturation_liquid"))->owner()) = 
	     *S_inter_->GetFieldData(getKey(name_ss.str(),"saturation_liquid"));

	   *S_next_->GetFieldData(getKey(name.str(),"unfrozen_fraction"), 
				   S_next_->GetField(getKey(name.str(),"unfrozen_fraction"))->owner()) = 
	     *S_inter_->GetFieldData(getKey(name.str(),"unfrozen_fraction"));
	   
	   *S_next_->GetFieldData(getKey(name.str(),"porosity"), 
				   S_next_->GetField(getKey(name.str(),"porosity"))->owner()) = 
	    *S_inter_->GetFieldData(getKey(name.str(),"porosity"));

	   *S_next_->GetFieldData(getKey(name.str(),"evaporative_flux"), 
				  S_next_->GetField(getKey(name.str(),"evaporative_flux"))->owner()) = 
	     *S_inter_->GetFieldData(getKey(name.str(),"evaporative_flux"));
	   
	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe1 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S_next_->GetFieldEvaluator(getKey(name.str(),"mass_source_temperature")));
	   
	   pfe1->SetFieldAsChanged(S_next_.ptr());

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe2 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S_next_->GetFieldEvaluator(getKey(name.str(),"conducted_energy_source")));
	   
	   pfe2->SetFieldAsChanged(S_next_.ptr());

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe3 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S_next_->GetFieldEvaluator(getKey(name.str(),"mass_source")));
	   
	   pfe3->SetFieldAsChanged(S_next_.ptr());

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe4 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S_next_->GetFieldEvaluator(getKey(name_ss.str(),"mass_source")));
	   
	   pfe4->SetFieldAsChanged(S_next_.ptr());


	   Teuchos::RCP<PK_PhysicalBDF_Default> pk_domain1 =
	     Teuchos::rcp_dynamic_cast<PK_PhysicalBDF_Default>(sub_pks_[count+1]);
	   ASSERT(pk_domain1.get()); // make sure the pk_domain is not empty
	   pk_domain1->ChangedSolution(S_next_.ptr());
	   
	   loc_dt = (*pk)->get_dt();
	   S_inter_->set_time(t0+t);
	   // revisit S_next_->set_time(t0 + t + loc_dt);
	}

      }
      count++;	
    }
    
   
  }
  

  MPI_Barrier(MPI_COMM_WORLD);  

  int nfailed_local = nfailed;
  S_->GetMesh("surface")->get_comm()->SumAll(&nfailed_local, &nfailed, 1);
 
 
  if (nfailed ==0){ 
    Epetra_MultiVector& surfstar_p = *S_next_->GetFieldData("surface_star-pressure",
							    S_inter_->GetField("surface_star-pressure")->owner())->ViewComponent("cell", false);
    Epetra_MultiVector& surfstar_t = *S_next_->GetFieldData("surface_star-temperature",
							    S_inter_->GetField("surface_star-temperature")->owner())->ViewComponent("cell", false);
    
    for (unsigned c=0; c<size_t; c++){
      std::stringstream name;
      int id = S_->GetMesh("surface")->cell_map(false).GID(c);
      name << "column_" << id <<"_surface";
      const Epetra_MultiVector& surf_p = *S_next_->GetFieldData(getKey(name.str(),"pressure"))->ViewComponent("cell", false);
      if(surf_p[0][0] > 101325.00)
	surfstar_p[0][c] = surf_p[0][0];
      else 
	surfstar_p[0][c]=101325.00;
      
    }
    
    for (unsigned c=0; c<size_t; c++){
      std::stringstream name;
      int id = S_->GetMesh("surface")->cell_map(false).GID(c);
      name << "column_" << id <<"_surface";
      const Epetra_MultiVector& surf_t = *S_next_->GetFieldData(getKey(name.str(),"temperature"))->ViewComponent("cell", false);
      surfstar_t[0][c] = surf_t[0][0];
    }

  
  // Mark surface_star-pressure evaluator as changed.
  // NOTE: later do it in the setup --aj
  Teuchos::RCP<PK_PhysicalBDF_Default> pk_surf =
    Teuchos::rcp_dynamic_cast<PK_PhysicalBDF_Default>(sub_pks_[0]);
  ASSERT(pk_surf.get());
  pk_surf->ChangedSolution();
  MPC<PK>::SubPKList::iterator pk1 = sub_pks_.begin();

   if(subcycle_key_)
     (*pk1)->CommitStep(t_old, t_new,S_next_);

  }
  if (nfailed > 0){
    flag_star = 1;
    return true;
  }
  else{
    flag_star = 1;
    return false;
  }
}
  

  bool WeakMPCSemiCoupled::CoupledSurfSubsurf3D(double t_old, double t_new, bool reinit) {
  bool fail = false;
  MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();

  // advance surface_star-pressure from t_n to t_(n+1)
  fail = (*pk)->AdvanceStep(t_old, t_new, reinit);
  
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
  
  Teuchos::RCP<PK_PhysicalBDF_Default> pk_domain =
    Teuchos::rcp_dynamic_cast<PK_PhysicalBDF_Default>(sub_pks_[1]);
  ASSERT(pk_domain.get()); // make sure the pk_domain is not empty
  pk_domain->ChangedSolution(S_inter_.ptr());
  
  if(fail) return fail;  
  // advance surface-pressure from t_n to t_(n+1)
  ++pk;
  fail += (*pk)->AdvanceStep(t_old, t_new, reinit);
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
  Teuchos::RCP<PK_PhysicalBDF_Default> pk_surf =
    Teuchos::rcp_dynamic_cast<PK_PhysicalBDF_Default>(sub_pks_[0]);
  ASSERT(pk_surf.get());
  pk_surf->ChangedSolution();
  
  return fail;
};



void 
WeakMPCSemiCoupled::generalize_inputspec(){
  
  Teuchos::Array<std::string> pks_order = plist_->get<Teuchos::Array<std::string> >("PKs order");
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_loc = S_loc->GetMesh("surface");
  int ncell = mesh_loc->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int pk_start = mesh_loc->cell_map(false).GID(0);
  int pk_end = mesh_loc->cell_map(false).GID(ncell-1) + 1;
  
  Teuchos::Array<std::string> pks_list_loc;
  pks_list_loc.push_back(pks_order[0]);
  pks_list_loc.push_back(pks_order[1]);
  for (int i = pk_start; i <pk_end; i++){
    std::stringstream name;
    name << "PK" << i;
    pks_list_loc.push_back(name.str());
    
  }
  plist_->set("PKs order", pks_list_loc);
  
  
  Teuchos::Array<std::string> pk_order = plist_->get<Teuchos::Array<std::string> >("PKs order"); // sublist of PK that correspond to each processor
  
  Teuchos::ParameterList pks_list_main = plist_->sublist("PKs");
  
    
  Teuchos::Array<std::string> pk2_order = pks_list_main.sublist(pk_order[1]).get<Teuchos::Array<std::string> >("PKs order"); // SEB PSS
  
    
  int len = pk_order.length();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
 // Generalizing PKs
  for(int i=2; i<len;i++){
    Teuchos::ParameterList pks_list = pks_list_main.sublist(pk_order[1]);
   
    pks_list.setName(pk_order[i]);
    
    std::string names_seb = pk2_order[0];
    Teuchos::Array<std::string> names = pks_list.sublist("PKs")
      .sublist(pk2_order[1]).get<Teuchos::Array<std::string> >("PKs order");

    std::string d_num = pk_order[i].substr(2,pk_order[i].size()-1);
    
    std::string domain_surf= "column_" + d_num + "_surface";
    std::string domain_ss = "column_"+ d_num;
    

    pks_list.sublist("PKs").sublist(names_seb)
      .set("primary variable",getKey(domain_surf,"snow_depth"));
    pks_list.sublist("PKs").sublist(names_seb)
      .set("conserved quantity key",getKey(domain_surf,"snow_depth"));
    pks_list.sublist("PKs").sublist(names_seb)
      .set("domain name",domain_surf);
    

    //PSS parameters
    
    pks_list.sublist("PKs").sublist(pk2_order[1])
      .sublist("ewc delegate").set("domain name",domain_ss);
    pks_list.sublist("PKs").sublist(pk2_order[1])
      .sublist("surface ewc delegate").set("domain name",domain_surf);
   
    Teuchos::ParameterList& pk2_list = pks_list.sublist("PKs").sublist(pk2_order[1]).sublist("PKs");

    pk2_list.sublist(names[0])
      .set("primary variable",getKey(domain_ss, "pressure"));
    pk2_list.sublist(names[0]).set("domain name",domain_ss);
    
    // if restarting from local checkpoint files
    
    if (pk2_list.sublist(names[0]).sublist("initial condition").isParameter("restart files, checkpoint cycles")){
   
      Teuchos::Array<std::string> restart = pk2_list.sublist(names[0]).sublist("initial condition").get<Teuchos::Array<std::string> >("restart files, checkpoint cycles");

      std::stringstream res_file;

      if(restart[0].rfind("/") == restart[0].size()-1){}
      else
	restart[0] += "/";
      
      res_file << restart[0] << "checkpoint_" << rank << "_" << restart[1] << ".h5";
    
      pk2_list.sublist(names[0]).sublist("initial condition").set("restart file", res_file.str());
    
    }
    
   
    
    pk2_list.sublist(names[0]).sublist("water retention evaluator")
      .set("rel perm key",getKey(domain_ss,"relative_permeability"));
    pk2_list.sublist(names[0]).sublist("water retention evaluator")
      .set("liquid saturation key",getKey(domain_ss,"saturation_liquid"));
    pk2_list.sublist(names[0]).sublist("water retention evaluator")
      .set("ice saturation key",getKey(domain_ss,"saturation_ice"));
    pk2_list.sublist(names[0]).sublist("water retention evaluator")
      .set("gas saturation key",getKey(domain_ss,"saturation_gas"));
    pk2_list.sublist(names[0]).sublist("water retention evaluator")
      .set("surface rel perm key",getKey(domain_surf,"relative_permeability"));
    
    
    pk2_list.sublist(names[1])
      .set("primary variable",getKey(domain_ss, "temperature"));
    pk2_list.sublist(names[1]).set("domain name",domain_ss);
    pk2_list.sublist(names[1]).sublist("thermal conductivity evaluator")
      .set("thermal conductivity key",getKey(domain_ss,"thermal_conductivity"));
    

    
      if (pk2_list.sublist(names[1]).sublist("initial condition").isParameter("restart files, checkpoint cycles")){
   
      Teuchos::Array<std::string> restart = pk2_list.sublist(names[1]).sublist("initial condition").get<Teuchos::Array<std::string> >("restart files, checkpoint cycles");

      std::stringstream res_file;

      if(restart[0].rfind("/") == restart[0].size()-1){}
      else
	restart[0] += "/";
      
      res_file << restart[0] << "checkpoint_" << rank << "_" << restart[1] << ".h5";
    
      pk2_list.sublist(names[1]).sublist("initial condition").set("restart file", res_file.str());
    
    }
    


    pk2_list.sublist(names[2])
      .set("primary variable",getKey(domain_surf, "pressure"));
    pk2_list.sublist(names[2])
      .set("domain name",domain_surf);
    pk2_list.sublist(names[2])
      .sublist("elevation evaluator").set("elevation key",getKey(domain_surf,"elevation"));
    pk2_list.sublist(names[2])
      .sublist("potential evaluator").set("potential key",getKey(domain_surf,"pres_elev"));
    pk2_list.sublist(names[2])
      .sublist("overland water content evaluator").set("domain name",domain_surf);
    pk2_list.sublist(names[2]).sublist("overland conductivity evaluator")
      .set("overland conductivity key",getKey(domain_surf,"overland_conductivity"));

    pk2_list.sublist(names[3])
      .set("primary variable",getKey(domain_surf, "temperature"));
    pk2_list.sublist(names[3])
      .set("domain name",domain_surf);
    pk2_list.sublist(names[3])
      .set("flux key",getKey(domain_surf, "mass_flux")); 
    pk2_list.sublist(names[3]).sublist("thermal conductivity evaluator")
      .set("thermal conductivity key",getKey(domain_surf,"thermal_conductivity"));


    plist_->sublist("PKs").set(pk_order[i], pks_list);
    
  }


  pk_order.remove(1); // we are done generalizing PKs (making multiple copies of 2nd PK in the inputspec with different domain names.) Lets delete it now
  plist_->set("PKs order", pk_order);


  Teuchos::ParameterList state_list =  FElist_loc;
  
  //Generalize state
  for(int i=1; i<len-1;i++){

    std::string d_num = pk_order[i].substr(2,pk_order[i].size()-1);
    std::string domain_surf= "column_" + d_num + "_surface";
    std::string domain_ss = "column_"+ d_num;
        
    Teuchos::ParameterList surf_wc = state_list.sublist("surface-water_content");
    surf_wc.setName(getKey(domain_surf,"water_content"));
    FElist_loc.set(surf_wc.name(), surf_wc);
 

    Teuchos::ParameterList surf_energy = state_list.sublist("surface-energy");
    surf_energy.setName(getKey(domain_surf,"energy"));
    FElist_loc.set(surf_energy.name(), surf_energy);
    

    Teuchos::ParameterList surf_pd = state_list.sublist("surface-ponded_depth");
    surf_pd.setName(getKey(domain_surf,"ponded_depth"));
    FElist_loc.set(surf_pd.name(), surf_pd);
 
    Teuchos::ParameterList surf_tes = state_list.sublist("surface-total_energy_source");
    surf_tes.setName(getKey(domain_surf,"total_energy_source"));
    surf_tes.set("domain",domain_surf);
    surf_tes.set("internal enthalpy key",getKey(domain_surf,"enthalpy"));
    surf_tes.set("external enthalpy key",getKey(domain_surf,"mass_source_enthalpy"));
    surf_tes.set("internal density key",getKey(domain_surf,"molar_density_liquid"));
    surf_tes.set("external density key",getKey(domain_surf,"source_molar_density"));
    FElist_loc.set(surf_tes.name(), surf_tes);

    Teuchos::ParameterList surf_mse = state_list.sublist("surface-mass_source_enthalpy");
    surf_mse.setName(getKey(domain_surf,"mass_source_enthalpy"));
    surf_mse.set("pressure key",getKey(domain_surf,"effective_pressure"));
    surf_mse.set("molar density key",getKey(domain_surf,"source_molar_density"));
    surf_mse.set("internal energy key",getKey(domain_surf,"source_internal_energy"));
    FElist_loc.set(surf_mse.name(), surf_mse);

    Teuchos::ParameterList surf_smd = state_list.sublist("surface-source_molar_density");
    surf_smd.setName(getKey(domain_surf,"source_molar_density"));
    surf_smd.set("pressure key",getKey(domain_surf,"effective_pressure"));
    surf_smd.set("molar density key",getKey(domain_surf,"source_molar_density"));
    surf_smd.set("temperature key",getKey(domain_surf,"mass_source_temperature"));
    FElist_loc.set(surf_smd.name(), surf_smd);

    Teuchos::ParameterList surf_sie = state_list.sublist("surface-source_internal_energy");
    surf_sie.setName(getKey(domain_surf,"source_internal_energy"));
    surf_sie.set("internal energy key",getKey(domain_surf,"source_internal_energy"));
    surf_sie.set("temperature key",getKey(domain_surf,"mass_source_temperature"));
    FElist_loc.set(surf_sie.name(), surf_sie);

    Teuchos::ParameterList surf_uf = state_list.sublist("surface-unfrozen_fraction");
    surf_uf.setName(getKey(domain_surf,"unfrozen_fraction"));
    FElist_loc.set(surf_uf.name(), surf_uf);

    Teuchos::ParameterList surf_mdl = state_list.sublist("surface-molar_density_liquid");
    surf_mdl.setName(getKey(domain_surf,"molar_density_liquid"));
    surf_mdl.set("pressure key",getKey(domain_surf,"effective_pressure"));
    surf_mdl.set("molar density key",getKey(domain_surf,"molar_density_liquid"));
    surf_mdl.set("mass density key",getKey(domain_surf,"mass_density_liquid"));
    FElist_loc.set(surf_mdl.name(), surf_mdl);

    Teuchos::ParameterList surf_ued = state_list.sublist("surface-unfrozen_effective_depth");
    surf_ued.setName(getKey(domain_surf,"unfrozen_effective_depth"));
    FElist_loc.set(surf_ued.name(), surf_ued);

    Teuchos::ParameterList surf_rp = state_list.sublist("surface-relative_permeability");
    surf_rp.setName(getKey(domain_surf,"relative_permeability"));
    FElist_loc.set(surf_rp.name(), surf_rp);


    Teuchos::ParameterList surf_mdi = state_list.sublist("surface-mass_density_ice");
    surf_mdi.setName(getKey(domain_surf,"mass_density_ice"));
    surf_mdi.set("pressure key",getKey(domain_surf,"effective_pressure"));
    surf_mdi.set("molar density key",getKey(domain_surf,"molar_density_ice"));
    surf_mdi.set("mass density key",getKey(domain_surf,"mass_density_ice"));
    FElist_loc.set(surf_mdi.name(), surf_mdi);

    Teuchos::ParameterList surf_modi = state_list.sublist("surface-molar_density_ice");
    surf_modi.setName(getKey(domain_surf,"molar_density_ice"));
    surf_modi.set("pressure key",getKey(domain_surf,"effective_pressure"));
    surf_modi.set("molar density key",getKey(domain_surf,"molar_density_ice"));
    surf_modi.set("mass density key",getKey(domain_surf,"mass_density_ice"));
    FElist_loc.set(surf_modi.name(), surf_modi);

    Teuchos::ParameterList surf_madl = state_list.sublist("surface-mass_density_liquid");
    surf_madl.setName(getKey(domain_surf,"mass_density_liquid"));
    surf_madl.set("pressure key",getKey(domain_surf,"effective_pressure"));
    surf_madl.set("molar density key",getKey(domain_surf,"molar_density_liquid"));
    surf_madl.set("mass density key",getKey(domain_surf,"mass_density_liquid"));
    FElist_loc.set(surf_madl.name(), surf_madl);

    Teuchos::ParameterList surf_iel = state_list.sublist("surface-internal_energy_liquid");
    surf_iel.setName(getKey(domain_surf,"internal_energy_liquid"));
    surf_iel.set("internal energy key",getKey(domain_surf,"internal_energy_liquid"));
    FElist_loc.set(surf_iel.name(), surf_iel);
 
    Teuchos::ParameterList surf_iei = state_list.sublist("surface-internal_energy_ice");
    surf_iei.setName(getKey(domain_surf,"internal_energy_ice"));
    surf_iei.set("internal energy key",getKey(domain_surf,"internal_energy_ice"));
    FElist_loc.set(surf_iei.name(), surf_iei);

    Teuchos::ParameterList surf_vp = state_list.sublist("surface-vapor_pressure");
    surf_vp.setName(getKey(domain_surf,"vapor_pressure"));
    surf_vp.set("surface key",getKey(domain_surf,"vapor_pressure"));
    surf_vp.set("subsurface key",getKey(domain_ss,"mol_frac_gas"));
    FElist_loc.set(surf_vp.name(), surf_vp);

    Teuchos::ParameterList surf_mc = state_list.sublist("surface-manning_coefficient");
    surf_mc.setName(getKey(domain_surf,"manning_coefficient"));
    FElist_loc.set(surf_mc.name(), surf_mc);

    Teuchos::ParameterList surf_ep = state_list.sublist("surface-effective_pressure");
    surf_ep.setName(getKey(domain_surf,"effective_pressure"));
    FElist_loc.set(surf_ep.name(), surf_ep);

    Teuchos::ParameterList surf_ilwr = state_list.sublist("surface-incoming_longwave_radiation");
    surf_ilwr.setName(getKey(domain_surf,"incoming_longwave_radiation"));
    FElist_loc.set(surf_ilwr.name(), surf_ilwr);
 

    Teuchos::ParameterList surf_iswr = state_list.sublist("surface-incoming_shortwave_radiation");
    surf_iswr.setName(getKey(domain_surf,"incoming_shortwave_radiation"));
    FElist_loc.set(surf_iswr.name(), surf_iswr);

    Teuchos::ParameterList surf_at = state_list.sublist("surface-air_temperature");
    surf_at.setName(getKey(domain_surf,"air_temperature"));
    FElist_loc.set(surf_at.name(), surf_at);

    Teuchos::ParameterList surf_rh = state_list.sublist("surface-relative_humidity");
    surf_rh.setName(getKey(domain_surf,"relative_humidity"));
    FElist_loc.set(surf_rh.name(), surf_rh);


    Teuchos::ParameterList surf_ws = state_list.sublist("surface-wind_speed");
    surf_ws.setName(getKey(domain_surf,"wind_speed"));
    FElist_loc.set(surf_ws.name(), surf_ws);
 
    Teuchos::ParameterList surf_pr = state_list.sublist("surface-precipitation_rain");
    surf_pr.setName(getKey(domain_surf,"precipitation_rain"));
    FElist_loc.set(surf_pr.name(), surf_pr);

    Teuchos::ParameterList surf_ps = state_list.sublist("surface-precipitation_snow");
    surf_ps.setName(getKey(domain_surf,"precipitation_snow"));
    FElist_loc.set(surf_ps.name(), surf_ps);

    Teuchos::ParameterList surf_sp = state_list.sublist("surface-porosity");
    surf_sp.setName(getKey(domain_surf,"porosity"));
    surf_sp.set("surface key",getKey(domain_surf,"porosity"));
    surf_sp.set("subsurface key",getKey(domain_ss,"porosity"));
    FElist_loc.set(surf_sp.name(), surf_sp);


    //---------------------- SUBSURFACE -----------
    
    Teuchos::ParameterList ss_wc = state_list.sublist("water_content");
    ss_wc.setName(getKey(domain_ss,"water_content"));
    FElist_loc.set(ss_wc.name(), ss_wc);

    Teuchos::ParameterList ss_energy = state_list.sublist("energy");
    ss_energy.setName(getKey(domain_ss,"energy"));
    FElist_loc.set(ss_energy.name(), ss_energy);
    
    Teuchos::ParameterList ss_cpgl = state_list.sublist("capillary_pressure_gas_liq");
    ss_cpgl.setName(getKey(domain_ss,"capillary_pressure_gas_liq"));
    FElist_loc.set(ss_cpgl.name(), ss_cpgl);

    Teuchos::ParameterList ss_cpil = state_list.sublist("capillary_pressure_liq_ice");
    ss_cpil.setName(getKey(domain_ss,"capillary_pressure_liq_ice"));
    FElist_loc.set(ss_cpil.name(), ss_cpil);

    Teuchos::ParameterList ss_mdl = state_list.sublist("molar_density_liquid");
    ss_mdl.setName(getKey(domain_ss,"molar_density_liquid"));
    ss_mdl.set("pressure key",getKey(domain_ss,"effective_pressure"));
    ss_mdl.set("molar density key",getKey(domain_ss,"molar_density_liquid"));
    ss_mdl.set("mass density key",getKey(domain_ss,"mass_density_liquid"));
    FElist_loc.set(ss_mdl.name(), ss_mdl);

    Teuchos::ParameterList ss_vis = state_list.sublist("viscosity_liquid");
    ss_vis.setName(getKey(domain_ss,"viscosity_liquid"));
    ss_vis.set("viscosity key",getKey(domain_ss,"viscosity_liquid"));
    FElist_loc.set(ss_vis.name(), ss_vis);

    Teuchos::ParameterList ss_mdg = state_list.sublist("molar_density_gas");
    ss_mdg.setName(getKey(domain_ss,"molar_density_gas"));
    ss_mdg.set("molar density key",getKey(domain_ss,"molar_density_gas"));
    FElist_loc.set(ss_mdg.name(), ss_mdg);

     
    Teuchos::ParameterList ss_mdi = state_list.sublist("molar_density_ice");
    ss_mdi.setName(getKey(domain_ss,"molar_density_ice"));
    ss_mdi.set("molar density key",getKey(domain_ss,"molar_density_ice"));
    FElist_loc.set(ss_mdi.name(), ss_mdi);


    Teuchos::ParameterList ss_iel = state_list.sublist("internal_energy_liquid");
    ss_iel.setName(getKey(domain_ss,"internal_energy_liquid"));
    ss_iel.set("molar density key",getKey(domain_ss,"internal_energy_liquid"));
    FElist_loc.set(ss_iel.name(), ss_iel);

    Teuchos::ParameterList ss_ier = state_list.sublist("internal_energy_rock");
    ss_ier.setName(getKey(domain_ss,"internal_energy_rock"));
    ss_ier.set("internal energy key",getKey(domain_ss,"internal_energy_rock"));
    FElist_loc.set(ss_ier.name(), ss_ier);

    Teuchos::ParameterList ss_ieg = state_list.sublist("internal_energy_gas");
    ss_ieg.setName(getKey(domain_ss,"internal_energy_gas"));
    ss_ieg.set("internal energy key",getKey(domain_ss,"internal_energy_gas"));
    FElist_loc.set(ss_ieg.name(), ss_ieg);

    Teuchos::ParameterList ss_iei = state_list.sublist("internal_energy_ice");
    ss_iei.setName(getKey(domain_ss,"internal_energy_ice"));
    ss_iei.set("internal energy key",getKey(domain_ss,"internal_energy_ice"));
    FElist_loc.set(ss_iei.name(), ss_iei);

    Teuchos::ParameterList ss_mf = state_list.sublist("mol_frac_gas");
    ss_mf.setName(getKey(domain_ss,"mol_frac_gas"));
    ss_mf.set("molar fraction key",getKey(domain_ss,"mol_frac_gas"));
    FElist_loc.set(ss_mf.name(), ss_mf);

    Teuchos::ParameterList ss_bp = state_list.sublist("base_porosity");
    ss_bp.setName(getKey(domain_ss,"base_porosity"));
    FElist_loc.set(ss_bp.name(), ss_bp);

    Teuchos::ParameterList ss_por = state_list.sublist("porosity");
    ss_por.setName(getKey(domain_ss,"porosity"));
    FElist_loc.set(ss_por.name(), ss_por);
    
    Teuchos::ParameterList ss_perm = state_list.sublist("permeability");
    ss_perm.setName(getKey(domain_ss,"permeability"));
    FElist_loc.set(ss_perm.name(), ss_perm);

    Teuchos::ParameterList ss_dr = state_list.sublist("density_rock");
    ss_dr.setName(getKey(domain_ss,"density_rock"));
    FElist_loc.set(ss_dr.name(), ss_dr);

    Teuchos::ParameterList ss_ep = state_list.sublist("effective_pressure");
    ss_ep.setName(getKey(domain_ss,"effective_pressure"));
    FElist_loc.set(ss_ep.name(), ss_ep);
    
  }
 
  
}
  
} // namespace Amanzi


