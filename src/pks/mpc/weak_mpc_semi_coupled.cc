#include "Teuchos_XMLParameterListHelpers.hpp"

//#include "pk_physical_bdf_base.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "strong_mpc.hh"

#include "weak_mpc_semi_coupled.hh"
#include "weak_mpc_semi_coupled_helper.hh"



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
  subcycle_key_ = plist_->get<bool>("subcycle",false);

  //Teuchos::ParameterList en_list = S->FEList().isSublist("surface_star-depression_depth");
  if (S->FEList().isSublist("surface_star-depression_depth"))
    sg_model_ = true;


  //  sg_model_ = en_list.get<bool>("subgrid model", false);
  if (sg_model_){
    delta_max_ = .4;//en_list.get<double>("maximum ponded depth");
    delta_ex_ = 0.2;//en_list.get<double>("excluded volume");
  }
  if(sg_model_)
    assert(delta_max_ > 0.0);
  //  sync_time_ = plist_->get<double>("sync time"); //provide default value later!!
  ASSERT(!(coupling_key_.empty()));

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
bool 
WeakMPCSemiCoupled::AdvanceStep(double t_old, double t_new, bool reinit) {
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
  double M_ = 0.0180153; //molar mass

  MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
  
  // advance surface_star-pressure from t_n to t_(n+1)

  //ensure the star solution is marked as changed when the subsurface columns fail
  if (flag_star || flag_star_surf){
    
    flag_star = 0;
    flag_star_surf=0;
    
    Teuchos::RCP<PK_BDF_Default> pk_sfstar =
      Teuchos::rcp_dynamic_cast<PK_BDF_Default>(sub_pks_[0]);
    
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

 
  const Epetra_MultiVector& surfstar_temp = *S_next_->GetFieldData("surface_star-temperature")
    ->ViewComponent("cell", false);
  unsigned int size_t = surfstar_temp.MyLength();

  assert(size_t == numPKs_ -1); // check if the subsurface columns are equal to the surface cells
  
  
  //copying pressure
  if(!sg_model_){
    const Epetra_MultiVector& surfstar_pres = *S_next_->GetFieldData("surface_star-pressure")->ViewComponent("cell", false);
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
  }
  else{
    
    const Epetra_MultiVector& vol_pd = *S_next_->GetFieldData("surface_star-volumetric_ponded_depth")
      ->ViewComponent("cell", false);
    
    const Epetra_MultiVector& mdl = *S_next_->GetFieldData("surface_star-mass_density_liquid")
      ->ViewComponent("cell", false);

    //const Epetra_MultiVector& surfstar_wc = *S_next_->GetFieldData("surface_star-water_content")
    // ->ViewComponent("cell", false);
    
    const Epetra_Vector& gravity = *S_->GetConstantVectorData("gravity");
    double gz = -gravity[2];
    
    for (unsigned c=0; c<size_t; c++){
    
      double pres = vol_pd[0][c]*mdl[0][c]*gz + 101325.0; // convert volumetric head to pressure
    
      if(pres > 101325.0){
	std::stringstream name;
	int id = S_->GetMesh("surface")->cell_map(false).GID(c);
	
	name << "column_" << id <<"_surface";
	Epetra_MultiVector& surf_pres = *S_inter_->GetFieldData(getKey(name.str(),"pressure"), 
								S_inter_->GetField(getKey(name.str(),"pressure"))->owner())
	  ->ViewComponent("cell", false);
	surf_pres[0][0] = pres;
      }
      else {}
    }
    
    
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
    Teuchos::RCP<PK_BDF_Default> pk_domain =
      Teuchos::rcp_dynamic_cast<PK_BDF_Default>(sub_pks_[i]);
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

	    UpdateIntermediateStateParameters(S_next_, S_inter_,id);

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe1 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
	     (S_inter_->GetFieldEvaluator(getKey(name.str(),"mass_source_temperature")));
	   
	   pfe1->SetFieldAsChanged(S_inter_.ptr());

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe2 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
	     (S_inter_->GetFieldEvaluator(getKey(name.str(),"conducted_energy_source")));
	   
	   pfe2->SetFieldAsChanged(S_inter_.ptr());

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe3 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
	     (S_inter_->GetFieldEvaluator(getKey(name.str(),"mass_source")));
	   
	   pfe3->SetFieldAsChanged(S_inter_.ptr());

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe4 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
	     (S_inter_->GetFieldEvaluator(getKey(name_ss.str(),"mass_source")));
	   
	   pfe4->SetFieldAsChanged(S_inter_.ptr());
	  
	   Teuchos::RCP<PK_PhysicalBDF_Default> pk_domain =
	     Teuchos::rcp_dynamic_cast<PK_PhysicalBDF_Default>(sub_pks_[count+1]);
	   ASSERT(pk_domain.get()); // make sure the pk_domain is not empty
	   pk_domain->ChangedSolution(S_inter_.ptr());
	   
	   S_inter_->set_time(t0+t);
	   S_next_->set_time( t0 + t + loc_dt);
	   
	  }
	  
	}
	else{
	  cyc_flag = true;
	  
	  UpdateNextStateParameters(S_next_, S_inter_, id);

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe1 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
	     (S_next_->GetFieldEvaluator(getKey(name.str(),"mass_source_temperature")));
	   
	   pfe1->SetFieldAsChanged(S_next_.ptr());

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe2 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
	     (S_next_->GetFieldEvaluator(getKey(name.str(),"conducted_energy_source")));
	   
	   pfe2->SetFieldAsChanged(S_next_.ptr());

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe3 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
	     (S_next_->GetFieldEvaluator(getKey(name.str(),"mass_source")));
	   
	   pfe3->SetFieldAsChanged(S_next_.ptr());

	   Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe4 = 
	     Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
	     (S_next_->GetFieldEvaluator(getKey(name_ss.str(),"mass_source")));
	   
	   pfe4->SetFieldAsChanged(S_next_.ptr());


	   Teuchos::RCP<PK_PhysicalBDF_Default> pk_domain =
	     Teuchos::rcp_dynamic_cast<PK_PhysicalBDF_Default>(sub_pks_[count+1]);
	   ASSERT(pk_domain.get()); // make sure the pk_domain is not empty
	   pk_domain->ChangedSolution(S_next_.ptr());
	   
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
    Epetra_MultiVector& surfstar_wc = *S_next_->GetFieldData("surface_star-water_content",
							     S_inter_->GetField("surface_star-water_content")->owner())->ViewComponent("cell", false);
    if (!sg_model_){
      for (unsigned c=0; c<size_t; c++){
	std::stringstream name;
	int id = S_->GetMesh("surface")->cell_map(false).GID(c);
	name << "column_" << id <<"_surface";
	const Epetra_MultiVector& surf_p = *S_next_->GetFieldData(getKey(name.str(),"pressure"))->ViewComponent("cell", false);
	const Epetra_MultiVector& surf_wc = *S_next_->GetFieldData(getKey(name.str(),"water_content"))->ViewComponent("cell", false);
	if(surf_p[0][0] > 101325.00){
	  surfstar_p[0][c] = surf_p[0][0];
	  surfstar_wc[0][c] = surf_wc[0][0];
	}
	else 
	  surfstar_p[0][c]=101325.00;	
      }
    }
    else{
      
      const Epetra_MultiVector& delta_max_v = *S_next_->GetFieldData("surface_star-maximum_ponded_depth")->ViewComponent("cell", false);
      const Epetra_MultiVector& delta_ex_v = *S_next_->GetFieldData("surface_star-excluded_volume")->ViewComponent("cell", false);
      //	const Epetra_MultiVector& depr_depth_v = *S->GetFieldData(depr_depth_key_)->ViewComponent("cell", false);
      const Epetra_Vector& gravity = *S_->GetConstantVectorData("gravity");
      double gz = -gravity[2];
      const double& p_atm = *S_->GetScalarData("atmospheric_pressure");

      for (unsigned c=0; c<size_t; c++){
	std::stringstream name;
	int id = S_->GetMesh("surface")->cell_map(false).GID(c);
	name << "column_" << id <<"_surface";
	const Epetra_MultiVector& pd = *S_next_->GetFieldData(getKey(name.str(),"ponded_depth"))->ViewComponent("cell", false);
	const Epetra_MultiVector& surf_wc = *S_next_->GetFieldData(getKey(name.str(),"water_content"))->ViewComponent("cell", false);
	
	const Epetra_MultiVector& cv = *S_next_->GetFieldData(getKey(name.str(),"cell_volume"))->ViewComponent("cell", false);

	const Epetra_MultiVector& mdl = *S_next_->GetFieldData(getKey(name.str(),"mass_density_liquid"))->ViewComponent("cell", false);
	
	
	if (pd[0][0] >0){
	  double delta = FindVolumetricHead(pd[0][0], delta_max_v[0][c],delta_ex_v[0][c]);
	  
	  double pres = delta*mdl[0][0] *gz + p_atm;
	  surfstar_p[0][c] = pres; 

	  double vpd = 0;
	  if (delta <= delta_max_v[0][c]){
	    vpd = std::pow(delta,2)*(2*delta_max_v[0][c] - 3*delta_ex_v[0][c])/std::pow(delta_max_v[0][c],2) 
	      + std::pow(delta,3)*(2*delta_ex_v[0][c] - delta_max_v[0][c])/std::pow(delta_max_v[0][c],3); //later call get the volumetric head given ponded depth to fix this
	  }
	  else{
	    vpd = delta - delta_ex_v[0][c];
	  }

	  double vpd_pres = vpd *mdl[0][0] *gz + p_atm;
	  
	  surfstar_wc[0][c] = (vpd_pres - p_atm)/ (gz * M_);
 	  surfstar_wc[0][c] *= cv[0][0];
	}
	else 
	  surfstar_p[0][c]=101325.0;
      }
      
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
  Teuchos::RCP<PK_BDF_Default> pk_surf =
    Teuchos::rcp_dynamic_cast<PK_BDF_Default>(sub_pks_[0]);
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
  

bool 
WeakMPCSemiCoupled::CoupledSurfSubsurf3D(double t_old, double t_new, bool reinit) {
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
WeakMPCSemiCoupled::generalize_inputspec(const Teuchos::Ptr<State>& S){
  
    
 //  Teuchos::Array<std::string> pks_order = plist_->get<Teuchos::Array<std::string> >("PKs order");
  Teuchos::ParameterList& top_pk_list = plist_->sublist("PKs");
  Teuchos::ParameterList& top_pk_mpc = top_pk_list.sublist("Top level MPC");
  Teuchos::ParameterList& cycle_driver_list = plist_->sublist("cycle driver").sublist("PK tree").sublist("Top level MPC");

  Teuchos::Array<std::string> pks_order = top_pk_mpc.get<Teuchos::Array<std::string> >("PKs order");
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_loc = S->GetMesh("surface");
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
  //  plist_->set("PKs order", pks_list_loc);
  top_pk_mpc.set("PKs order", pks_list_loc);
  
 
  Teuchos::Array<std::string> pk_order = top_pk_mpc.get<Teuchos::Array<std::string> >("PKs order"); // sublist of PK that correspond to each processor
  
  Teuchos::ParameterList pks_list_main = top_pk_list;//plist_->sublist("PKs");
  
  
  Teuchos::Array<std::string> pk2_order = top_pk_list.sublist(pk_order[1]).get<Teuchos::Array<std::string> >("PKs order"); // SEB PSS

    
  int len = pk_order.length();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
 // Generalizing PKs

  for(int i=2; i<len;i++){
    std::string d_num = pk_order[i].substr(2,pk_order[i].size()-1);
    
    Teuchos::ParameterList temp;
    temp.setName(pk_order[i]);
    temp.set("PK type", "strong MPC");
    Teuchos::ParameterList& temp1 = temp.sublist("SEB"+d_num);
    temp1.set("PK type", "surface balance implicit");
    
    Teuchos::ParameterList& pss_list_cyc = temp.sublist("PSS"+d_num);
    pss_list_cyc.set("PK type", "permafrost model");

    Teuchos::ParameterList& temp2 = pss_list_cyc.sublist("SSF"+d_num);
    temp2.set("PK type", "permafrost flow");


    Teuchos::ParameterList& temp3 = pss_list_cyc.sublist("SSE"+d_num);
    temp3.set("PK type", "three-phase energy");
    
    Teuchos::ParameterList& temp4 = pss_list_cyc.sublist("SF"+d_num);
    temp4.set("PK type", "overland flow with ice");

    Teuchos::ParameterList& temp5 = pss_list_cyc.sublist("SE"+d_num);
    temp5.set("PK type", "surface energy");

    Teuchos::ParameterList sub_cyc = cycle_driver_list.sublist(pk_order[i]); 
    sub_cyc = temp;

    cycle_driver_list.set(pk_order[i],sub_cyc);

    FElist_loc.set(pk_order[i],temp);
   
    
   

    
    Teuchos::ParameterList pks_list = pks_list_main.sublist(pk_order[1]); //PKG2 list
   
    pks_list.setName(pk_order[i]); // change the name of PKG2 to column version
   
   
    pks_list.set("PK name",pk_order[i]);
    std::string domain_surf= "column_" + d_num + "_surface";
    std::string domain_ss = "column_"+ d_num;
    
    Teuchos::Array<std::string> seb_pss;
    std::string seb_pk = pk2_order[0] + d_num;
    std::string pss_pk = pk2_order[1] + d_num;
    seb_pss.push_back(seb_pk);
    seb_pss.push_back(pss_pk);
    pks_list.set("PKs order", seb_pss);

    top_pk_list.set(pk_order[i], pks_list);
    //    plist_->sublist("PKs").set(pk_order[i], pks_list);
   
      ///------
   

   
    Teuchos::ParameterList seb_list = pks_list_main.sublist(pk2_order[0]); //SEB list
    seb_list.setName(seb_pk); // change the name of SEB list
    std::string names_seb = pk2_order[0] +d_num;

    seb_list.set("primary variable",getKey(domain_surf,"snow_depth"));
    seb_list.set("conserved quantity key",getKey(domain_surf,"snow_depth"));
    seb_list.set("domain name",domain_surf);
    

    if (seb_list.sublist("initial condition").isParameter("restart files, checkpoint cycles")){
   
      Teuchos::Array<std::string> restart = seb_list.sublist("initial condition").get<Teuchos::Array<std::string> >("restart files, checkpoint cycles");

      std::stringstream res_file;

      if(restart[0].rfind("/") == restart[0].size()-1){}
      else
	restart[0] += "/";
      
      res_file << restart[0] << "checkpoint_" << rank << "_" << restart[1] << ".h5";
    
      seb_list.sublist("initial condition").set("restart file", res_file.str());
    
    }


    top_pk_list.set(seb_pk, seb_list);

    //------------------- PSS
    Teuchos::ParameterList pss_list = pks_list_main.sublist(pk2_order[1]); //PSS list
    pss_list.setName(pss_pk); // change the name of PSS lists

    pss_list.set("PK name",pss_pk);

    Teuchos::Array<std::string> names = pks_list_main
      .sublist(pk2_order[1]).get<Teuchos::Array<std::string> >("PKs order"); // get the PK order in PSS={SSF,SSE,SF,SE}
    Teuchos::Array<std::string> pks_pss;
    std::string ssf = names[0] + d_num;
    std::string sse = names[1] + d_num;
    std::string sf = names[2] + d_num;
    std::string se = names[3] + d_num;
    pks_pss.push_back(ssf);
    pks_pss.push_back(sse);
    pks_pss.push_back(sf);
    pks_pss.push_back(se);

    pss_list.set("PKs order", pks_pss);

    
    //PSS parameters
    
    pss_list.sublist("ewc delegate").set("domain name",domain_ss);
    pss_list.sublist("surface ewc delegate").set("domain name",domain_surf);
   
    
    top_pk_list.set(pss_pk, pss_list);

    // -------------- SSF
    Teuchos::ParameterList ssf_list = pks_list_main.sublist(names[0]);

    ssf_list.set("primary variable",getKey(domain_ss, "pressure"));
    ssf_list.set("domain name",domain_ss);
    
    // if restarting from local checkpoint files
    
    if (ssf_list.sublist("initial condition").isParameter("restart files, checkpoint cycles")){
   
      Teuchos::Array<std::string> restart = ssf_list.sublist("initial condition").get<Teuchos::Array<std::string> >("restart files, checkpoint cycles");

      std::stringstream res_file;

      if(restart[0].rfind("/") == restart[0].size()-1){}
      else
	restart[0] += "/";
      
      res_file << restart[0] << "checkpoint_" << rank << "_" << restart[1] << ".h5";
    
      ssf_list.sublist("initial condition").set("restart file", res_file.str());
    
    }
    
   
    
    ssf_list.sublist("water retention evaluator")
      .set("rel perm key",getKey(domain_ss,"relative_permeability"));
    ssf_list.sublist("water retention evaluator")
      .set("liquid saturation key",getKey(domain_ss,"saturation_liquid"));
    ssf_list.sublist("water retention evaluator")
      .set("ice saturation key",getKey(domain_ss,"saturation_ice"));
    ssf_list.sublist("water retention evaluator")
      .set("gas saturation key",getKey(domain_ss,"saturation_gas"));
    ssf_list.sublist("water retention evaluator")
      .set("surface rel perm key",getKey(domain_surf,"relative_permeability"));
    
    
    top_pk_list.set(ssf, ssf_list);
    

    // -------------- SSE
    Teuchos::ParameterList sse_list = pks_list_main.sublist(names[1]);
    
    sse_list.set("primary variable",getKey(domain_ss, "temperature"));
    sse_list.set("domain name",domain_ss);
    sse_list.sublist("thermal conductivity evaluator")
      .set("thermal conductivity key",getKey(domain_ss,"thermal_conductivity"));
    
    
    
    if (sse_list.sublist("initial condition").isParameter("restart files, checkpoint cycles")){
      
      Teuchos::Array<std::string> restart = sse_list.sublist("initial condition").get<Teuchos::Array<std::string> >("restart files, checkpoint cycles");
      
      std::stringstream res_file;
      
      if(restart[0].rfind("/") == restart[0].size()-1){}
      else
	restart[0] += "/";
      
      res_file << restart[0] << "checkpoint_" << rank << "_" << restart[1] << ".h5";
      
      sse_list.sublist("initial condition").set("restart file", res_file.str());
    
    }
    
    top_pk_list.set(sse, sse_list);
     // -------------- SSE
    Teuchos::ParameterList sf_list = pks_list_main.sublist(names[2]);

    sf_list.set("primary variable",getKey(domain_surf, "pressure"));
    sf_list.set("domain name",domain_surf);
    sf_list.sublist("elevation evaluator").set("elevation key",getKey(domain_surf,"elevation"));
    sf_list.sublist("potential evaluator").set("potential key",getKey(domain_surf,"pres_elev"));
    sf_list.sublist("overland water content evaluator").set("domain name",domain_surf);
    sf_list.sublist("overland conductivity evaluator")
      .set("overland conductivity key",getKey(domain_surf,"overland_conductivity"));

    // -------------- SE
    Teuchos::ParameterList se_list = pks_list_main.sublist(names[3]);
    se_list.set("primary variable",getKey(domain_surf, "temperature"));
    se_list.set("domain name",domain_surf);
    se_list.set("flux key",getKey(domain_surf, "mass_flux")); 
    se_list.sublist("thermal conductivity evaluator")
      .set("thermal conductivity key",getKey(domain_surf,"thermal_conductivity"));
    
    top_pk_list.set(sf, sf_list);
    top_pk_list.set(se, se_list);

    //    plist_->sublist("PKs").set(pk_order[i], pks_list);

  }


  pk_order.remove(1); // we are done generalizing PKs (making multiple copies of 2nd PK in the inputspec with different domain names.) Lets delete it now
  //  plist_->set("PKs order", pk_order);
  top_pk_mpc.set("PKs order", pk_order);

  
  Teuchos::ParameterList& state_list = S->FEList();//plist_->sublist("state").sublist("field evaluators");
 
  Teuchos::ParameterList state_list_loc = state_list;
  //Generalize state
  for(int i=1; i<len-1;i++){

    std::string d_num = pk_order[i].substr(2,pk_order[i].size()-1);
    std::string domain_surf= "column_" + d_num + "_surface";
    std::string domain_ss = "column_"+ d_num;
        
    Teuchos::ParameterList surf_wc = state_list_loc.sublist("surface-water_content");
    surf_wc.setName(getKey(domain_surf,"water_content"));
    state_list.set(surf_wc.name(), surf_wc);
    

    Teuchos::ParameterList surf_energy = state_list_loc.sublist("surface-energy");
    surf_energy.setName(getKey(domain_surf,"energy"));
    state_list.set(surf_energy.name(), surf_energy);
    
    
    Teuchos::ParameterList surf_pd = state_list_loc.sublist("surface-ponded_depth");
    surf_pd.setName(getKey(domain_surf,"ponded_depth"));
    state_list.set(surf_pd.name(), surf_pd);
 
    Teuchos::ParameterList surf_tes = state_list_loc.sublist("surface-total_energy_source");
    surf_tes.setName(getKey(domain_surf,"total_energy_source"));
    surf_tes.set("domain",domain_surf);
    surf_tes.set("internal enthalpy key",getKey(domain_surf,"enthalpy"));
    surf_tes.set("external enthalpy key",getKey(domain_surf,"mass_source_enthalpy"));
    surf_tes.set("internal density key",getKey(domain_surf,"molar_density_liquid"));
    surf_tes.set("external density key",getKey(domain_surf,"source_molar_density"));
    state_list.set(surf_tes.name(), surf_tes);

    Teuchos::ParameterList surf_mse = state_list_loc.sublist("surface-mass_source_enthalpy");
    surf_mse.setName(getKey(domain_surf,"mass_source_enthalpy"));
    surf_mse.set("pressure key",getKey(domain_surf,"effective_pressure"));
    surf_mse.set("molar density key",getKey(domain_surf,"source_molar_density"));
    surf_mse.set("internal energy key",getKey(domain_surf,"source_internal_energy"));
    state_list.set(surf_mse.name(), surf_mse);

    Teuchos::ParameterList surf_smd = state_list_loc.sublist("surface-source_molar_density");
    surf_smd.setName(getKey(domain_surf,"source_molar_density"));
    surf_smd.set("pressure key",getKey(domain_surf,"effective_pressure"));
    surf_smd.set("molar density key",getKey(domain_surf,"source_molar_density"));
    surf_smd.set("temperature key",getKey(domain_surf,"mass_source_temperature"));
    state_list.set(surf_smd.name(), surf_smd);

    Teuchos::ParameterList surf_sie = state_list_loc.sublist("surface-source_internal_energy");
    surf_sie.setName(getKey(domain_surf,"source_internal_energy"));
    surf_sie.set("internal energy key",getKey(domain_surf,"source_internal_energy"));
    surf_sie.set("temperature key",getKey(domain_surf,"mass_source_temperature"));
    state_list.set(surf_sie.name(), surf_sie);

    Teuchos::ParameterList surf_uf = state_list_loc.sublist("surface-unfrozen_fraction");
    surf_uf.setName(getKey(domain_surf,"unfrozen_fraction"));
    state_list.set(surf_uf.name(), surf_uf);

    Teuchos::ParameterList surf_mdl = state_list_loc.sublist("surface-molar_density_liquid");
    surf_mdl.setName(getKey(domain_surf,"molar_density_liquid"));
    surf_mdl.set("pressure key",getKey(domain_surf,"effective_pressure"));
    surf_mdl.set("molar density key",getKey(domain_surf,"molar_density_liquid"));
    surf_mdl.set("mass density key",getKey(domain_surf,"mass_density_liquid"));
    state_list.set(surf_mdl.name(), surf_mdl);

    Teuchos::ParameterList surf_ued = state_list_loc.sublist("surface-unfrozen_effective_depth");
    surf_ued.setName(getKey(domain_surf,"unfrozen_effective_depth"));
    state_list.set(surf_ued.name(), surf_ued);

    Teuchos::ParameterList surf_rp = state_list_loc.sublist("surface-relative_permeability");
    surf_rp.setName(getKey(domain_surf,"relative_permeability"));
    state_list.set(surf_rp.name(), surf_rp);


    Teuchos::ParameterList surf_mdi = state_list_loc.sublist("surface-mass_density_ice");
    surf_mdi.setName(getKey(domain_surf,"mass_density_ice"));
    surf_mdi.set("pressure key",getKey(domain_surf,"effective_pressure"));
    surf_mdi.set("molar density key",getKey(domain_surf,"molar_density_ice"));
    surf_mdi.set("mass density key",getKey(domain_surf,"mass_density_ice"));
    state_list.set(surf_mdi.name(), surf_mdi);

    Teuchos::ParameterList surf_modi = state_list_loc.sublist("surface-molar_density_ice");
    surf_modi.setName(getKey(domain_surf,"molar_density_ice"));
    surf_modi.set("pressure key",getKey(domain_surf,"effective_pressure"));
    surf_modi.set("molar density key",getKey(domain_surf,"molar_density_ice"));
    surf_modi.set("mass density key",getKey(domain_surf,"mass_density_ice"));
    state_list.set(surf_modi.name(), surf_modi);

    Teuchos::ParameterList surf_madl = state_list_loc.sublist("surface-mass_density_liquid");
    surf_madl.setName(getKey(domain_surf,"mass_density_liquid"));
    surf_madl.set("pressure key",getKey(domain_surf,"effective_pressure"));
    surf_madl.set("molar density key",getKey(domain_surf,"molar_density_liquid"));
    surf_madl.set("mass density key",getKey(domain_surf,"mass_density_liquid"));
    state_list.set(surf_madl.name(), surf_madl);

    Teuchos::ParameterList surf_iel = state_list_loc.sublist("surface-internal_energy_liquid");
    surf_iel.setName(getKey(domain_surf,"internal_energy_liquid"));
    surf_iel.set("internal energy key",getKey(domain_surf,"internal_energy_liquid"));
    state_list.set(surf_iel.name(), surf_iel);
 
    Teuchos::ParameterList surf_iei = state_list_loc.sublist("surface-internal_energy_ice");
    surf_iei.setName(getKey(domain_surf,"internal_energy_ice"));
    surf_iei.set("internal energy key",getKey(domain_surf,"internal_energy_ice"));
    state_list.set(surf_iei.name(), surf_iei);

    Teuchos::ParameterList surf_vp = state_list_loc.sublist("surface-vapor_pressure");
    surf_vp.setName(getKey(domain_surf,"vapor_pressure"));
    surf_vp.set("surface key",getKey(domain_surf,"vapor_pressure"));
    surf_vp.set("subsurface key",getKey(domain_ss,"mol_frac_gas"));
    state_list.set(surf_vp.name(), surf_vp);

    Teuchos::ParameterList surf_mc = state_list_loc.sublist("surface-manning_coefficient");
    surf_mc.setName(getKey(domain_surf,"manning_coefficient"));
    state_list.set(surf_mc.name(), surf_mc);

    Teuchos::ParameterList surf_ep = state_list_loc.sublist("surface-effective_pressure");
    surf_ep.setName(getKey(domain_surf,"effective_pressure"));
    state_list.set(surf_ep.name(), surf_ep);

    Teuchos::ParameterList surf_ilwr = state_list_loc.sublist("surface-incoming_longwave_radiation");
    surf_ilwr.setName(getKey(domain_surf,"incoming_longwave_radiation"));
    state_list.set(surf_ilwr.name(), surf_ilwr);
 

    Teuchos::ParameterList surf_iswr = state_list_loc.sublist("surface-incoming_shortwave_radiation");
    surf_iswr.setName(getKey(domain_surf,"incoming_shortwave_radiation"));
    state_list.set(surf_iswr.name(), surf_iswr);

    Teuchos::ParameterList surf_at = state_list_loc.sublist("surface-air_temperature");
    surf_at.setName(getKey(domain_surf,"air_temperature"));
    state_list.set(surf_at.name(), surf_at);

    Teuchos::ParameterList surf_rh = state_list_loc.sublist("surface-relative_humidity");
    surf_rh.setName(getKey(domain_surf,"relative_humidity"));
    state_list.set(surf_rh.name(), surf_rh);


    Teuchos::ParameterList surf_ws = state_list_loc.sublist("surface-wind_speed");
    surf_ws.setName(getKey(domain_surf,"wind_speed"));
    state_list.set(surf_ws.name(), surf_ws);
 
    Teuchos::ParameterList surf_pr = state_list_loc.sublist("surface-precipitation_rain");
    surf_pr.setName(getKey(domain_surf,"precipitation_rain"));
    state_list.set(surf_pr.name(), surf_pr);

    Teuchos::ParameterList surf_ps = state_list_loc.sublist("surface-precipitation_snow");
    surf_ps.setName(getKey(domain_surf,"precipitation_snow"));
    state_list.set(surf_ps.name(), surf_ps);

    Teuchos::ParameterList surf_sp = state_list_loc.sublist("surface-porosity");
    surf_sp.setName(getKey(domain_surf,"porosity"));
    surf_sp.set("surface key",getKey(domain_surf,"porosity"));
    surf_sp.set("subsurface key",getKey(domain_ss,"porosity"));
    state_list.set(surf_sp.name(), surf_sp);


    //---------------------- SUBSURFACE -----------
    
    Teuchos::ParameterList ss_wc = state_list_loc.sublist("water_content");
    ss_wc.setName(getKey(domain_ss,"water_content"));
    state_list.set(ss_wc.name(), ss_wc);

    Teuchos::ParameterList ss_energy = state_list_loc.sublist("energy");
    ss_energy.setName(getKey(domain_ss,"energy"));
    state_list.set(ss_energy.name(), ss_energy);
    
    Teuchos::ParameterList ss_cpgl = state_list_loc.sublist("capillary_pressure_gas_liq");
    ss_cpgl.setName(getKey(domain_ss,"capillary_pressure_gas_liq"));
    state_list.set(ss_cpgl.name(), ss_cpgl);

    Teuchos::ParameterList ss_cpil = state_list_loc.sublist("capillary_pressure_liq_ice");
    ss_cpil.setName(getKey(domain_ss,"capillary_pressure_liq_ice"));
    state_list.set(ss_cpil.name(), ss_cpil);

    Teuchos::ParameterList ss_mdl = state_list_loc.sublist("molar_density_liquid");
    ss_mdl.setName(getKey(domain_ss,"molar_density_liquid"));
    ss_mdl.set("pressure key",getKey(domain_ss,"effective_pressure"));
    ss_mdl.set("molar density key",getKey(domain_ss,"molar_density_liquid"));
    ss_mdl.set("mass density key",getKey(domain_ss,"mass_density_liquid"));
    state_list.set(ss_mdl.name(), ss_mdl);

    Teuchos::ParameterList ss_vis = state_list_loc.sublist("viscosity_liquid");
    ss_vis.setName(getKey(domain_ss,"viscosity_liquid"));
    ss_vis.set("viscosity key",getKey(domain_ss,"viscosity_liquid"));
    state_list.set(ss_vis.name(), ss_vis);

    Teuchos::ParameterList ss_mdg = state_list_loc.sublist("molar_density_gas");
    ss_mdg.setName(getKey(domain_ss,"molar_density_gas"));
    ss_mdg.set("molar density key",getKey(domain_ss,"molar_density_gas"));
    state_list.set(ss_mdg.name(), ss_mdg);

     
    Teuchos::ParameterList ss_mdi = state_list_loc.sublist("molar_density_ice");
    ss_mdi.setName(getKey(domain_ss,"molar_density_ice"));
    ss_mdi.set("molar density key",getKey(domain_ss,"molar_density_ice"));
    state_list.set(ss_mdi.name(), ss_mdi);


    Teuchos::ParameterList ss_iel = state_list_loc.sublist("internal_energy_liquid");
    ss_iel.setName(getKey(domain_ss,"internal_energy_liquid"));
    ss_iel.set("molar density key",getKey(domain_ss,"internal_energy_liquid"));
    state_list.set(ss_iel.name(), ss_iel);

    Teuchos::ParameterList ss_ier = state_list_loc.sublist("internal_energy_rock");
    ss_ier.setName(getKey(domain_ss,"internal_energy_rock"));
    ss_ier.set("internal energy key",getKey(domain_ss,"internal_energy_rock"));
    state_list.set(ss_ier.name(), ss_ier);

    Teuchos::ParameterList ss_ieg = state_list_loc.sublist("internal_energy_gas");
    ss_ieg.setName(getKey(domain_ss,"internal_energy_gas"));
    ss_ieg.set("internal energy key",getKey(domain_ss,"internal_energy_gas"));
    state_list.set(ss_ieg.name(), ss_ieg);

    Teuchos::ParameterList ss_iei = state_list_loc.sublist("internal_energy_ice");
    ss_iei.setName(getKey(domain_ss,"internal_energy_ice"));
    ss_iei.set("internal energy key",getKey(domain_ss,"internal_energy_ice"));
    state_list.set(ss_iei.name(), ss_iei);

    Teuchos::ParameterList ss_mf = state_list_loc.sublist("mol_frac_gas");
    ss_mf.setName(getKey(domain_ss,"mol_frac_gas"));
    ss_mf.set("molar fraction key",getKey(domain_ss,"mol_frac_gas"));
    state_list.set(ss_mf.name(), ss_mf);

    Teuchos::ParameterList ss_bp = state_list_loc.sublist("base_porosity");
    ss_bp.setName(getKey(domain_ss,"base_porosity"));
    state_list.set(ss_bp.name(), ss_bp);

    Teuchos::ParameterList ss_por = state_list_loc.sublist("porosity");
    ss_por.setName(getKey(domain_ss,"porosity"));
    state_list.set(ss_por.name(), ss_por);
    
    Teuchos::ParameterList ss_perm = state_list_loc.sublist("permeability");
    ss_perm.setName(getKey(domain_ss,"permeability"));
    state_list.set(ss_perm.name(), ss_perm);

    Teuchos::ParameterList ss_dr = state_list_loc.sublist("density_rock");
    ss_dr.setName(getKey(domain_ss,"density_rock"));
    state_list.set(ss_dr.name(), ss_dr);

    Teuchos::ParameterList ss_ep = state_list_loc.sublist("effective_pressure");
    ss_ep.setName(getKey(domain_ss,"effective_pressure"));
    state_list.set(ss_ep.name(), ss_ep);
    
  }
  
  plist_->sublist("state").set("field evaluators", state_list);
  
}




double 
WeakMPCSemiCoupled::FindVolumetricHead(double d, double delta_max, double delta_ex){

  double a = (2*delta_ex - delta_max) / std::pow(delta_max,3);
  double b = (2*delta_max - 3*delta_ex) / std::pow(delta_max,2);

  double x1=0,x2=delta_max,x3;
  int count=0;
  double tol = 1.0E-15;
  if (d <= delta_max){
    while (count <3000){
      x3 = (x1+x2)*0.5;
      double a1 = VolumetricHead(x1,a,b,d);
      double a3 = VolumetricHead(x3,a,b,d);
      
      if (a1*a3 <0)
	x2 = x3;
      else
	x1 = x3;
      if (std::fabs(VolumetricHead(x3,a,b,d)) < tol){
	break;
      }
      else
	count++;
    }
  }
  else
    return d + delta_ex;

  return x3;
}
double
WeakMPCSemiCoupled::VolumetricHead(double x, double a, double b, double d)
{
  double r = a*std::pow(x,3) + b*std::pow(x,2) - d;
  return r;
}
  
} // namespace Amanzi


