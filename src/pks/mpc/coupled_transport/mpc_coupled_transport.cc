/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

*/

#include "mpc_coupled_transport.hh"
#include "Transport_PK_ATS.hh"

namespace Amanzi {

  CoupledTransport_PK::CoupledTransport_PK(Teuchos::ParameterList& pk_tree_or_fe_list,
                                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                           const Teuchos::RCP<State>& S,
                                           const Teuchos::RCP<TreeVector>& soln) :
    PK_MPCSubcycled_ATS(pk_tree_or_fe_list, global_list, S, soln)
  {

    //std::cout << *global_list<<"\n";
    //exit(0);
      // Create verbosity object.
    vo_ = Teuchos::null;
    Teuchos::ParameterList vlist;
    vlist.sublist("verbose object") = global_list -> sublist("verbose object");
    vo_ =  Teuchos::rcp(new VerboseObject("Coupled TransportPK", vlist)); 

    Teuchos::OSTab tab = vo_->getOSTab();

  }

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double CoupledTransport_PK::get_dt() {
  //double dt = Amanzi::PK_MPCSubcycled_ATS::get_dt();
  master_dt_ = sub_pks_[master_]->get_dt();
  slave_dt_ = sub_pks_[slave_]->get_dt();

  std::cout<<"master dt = "<<master_dt_<<"\n";
  std::cout<<"slave_dt_ ="<<slave_dt_<<"\n";

  double dt = std::min(sub_pks_[master_]->get_dt(), sub_pks_[slave_]->get_dt());

  set_dt(dt);
  return dt;
}


// -----------------------------------------------------------------------------
// Set master dt
// -----------------------------------------------------------------------------
void CoupledTransport_PK::set_dt(double dt) {
  master_dt_ = dt;
  sub_pks_[master_]->set_dt(dt);
}


// -----------------------------------------------------------------------------
// Make necessary operatios by the end of the time steps.
// -----------------------------------------------------------------------------
void CoupledTransport_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {
  sub_pks_[master_]->CommitStep(t_old, t_new, S);
  sub_pks_[slave_]->CommitStep(t_old, t_new, S);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool CoupledTransport_PK::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;

  double dt_MPC = S_->final_time() - S_->initial_time();

  Teuchos::RCP<const CompositeVector> pond_prev = S_->GetFieldData("prev_ponded_depth");
  Teuchos::RCP<const CompositeVector> pond = S_next_->GetFieldData("ponded_depth");
  Teuchos::RCP<CompositeVector> pond_start  = S_->GetFieldCopyData("ponded_depth", "subcycle_start", "ponded_depth");
  if ((t_old > S_->initial_time())||(t_new < S_->final_time())) {
    InterpolateCellVector(*pond_prev->ViewComponent("cell"),
                          *pond->ViewComponent("cell"),
                          S_->initial_time() - t_old,
                          dt_MPC,
                          *pond_start->ViewComponent("cell"));
  }else{
    *pond_start = *pond_prev;   
    //S_->CopyField("ponded_depth", "subcycle_end", "ponded_depth");
  }

  Teuchos::RCP<const CompositeVector> sat_prev = S_->GetFieldData("prev_saturation_liquid");
  Teuchos::RCP<const CompositeVector> sat = S_next_->GetFieldData("saturation_liquid");
  Teuchos::RCP<CompositeVector> sat_start  = S_->GetFieldCopyData("saturation_liquid", "subcycle_start", "saturation_liquid");

  if ((t_old > S_->initial_time())||(t_new < S_->final_time())) {
    InterpolateCellVector(*sat_prev->ViewComponent("cell"),
                          *sat->ViewComponent("cell"),
                          S_->initial_time() - t_old,
                          dt_MPC,
                          *sat_start->ViewComponent("cell"));
  }else{
    *sat_start->ViewComponent("cell") = *sat_prev->ViewComponent("cell");
    //S_->CopyField("saturation_liquid","subcycle_end","saturation_liquid");
  }


  Teuchos::RCP<const CompositeVector> next_darcy = S_next_->GetFieldData("darcy_flux");
  Key flux_owner = S_next_->GetField("darcy_flux")->owner();
  Teuchos::RCP<CompositeVector>  next_darcy_copy = S_->GetFieldCopyData("darcy_flux", "next_timestep", flux_owner);
  *next_darcy_copy = *next_darcy;

  Teuchos::RCP<const CompositeVector> next_sur_flux = S_next_->GetFieldData("surface-flux");
  flux_owner = S_next_->GetField("surface-flux")->owner();
  Teuchos::RCP<CompositeVector>  next_sur_flux_copy = S_->GetFieldCopyData("surface-flux", "next_timestep", flux_owner);
  *next_sur_flux_copy = *next_sur_flux;

  sub_pks_[master_]->AdvanceStep(t_old, t_new, reinit);
  sub_pks_[slave_]->AdvanceStep(t_old, t_new, reinit);


  // const Epetra_MultiVector& surf_tcc = *S_->GetFieldCopyData("surface-total_component_concentration", "subcycling")->ViewComponent("cell",false);  
  // const Epetra_MultiVector& tcc = *S_->GetFieldCopyData("total_component_concentration", "subcycling")->ViewComponent("cell",false);

  // int num_components = 1;
  // std::vector<double> mass_subsurface(num_components, 0.), mass_surface(num_components, 0.);

  // for (int i=0; i<num_components; i++){
  //   mass_subsurface[i] = Teuchos::rcp_dynamic_cast<Transport::Transport_PK_ATS>(sub_pks_[master_])->ComputeSolute(tcc, i);
  //   mass_surface[i] = Teuchos::rcp_dynamic_cast<Transport::Transport_PK_ATS>(sub_pks_[slave_])->ComputeSolute(surf_tcc, i);

  //   *vo_->os() <<", subsurface =" << mass_subsurface[i] << " mol";
  //   *vo_->os() <<", surface =" << mass_surface[i]<< " mol";
  //   *vo_->os() <<", ToTaL =" << mass_surface[i]+mass_subsurface[i]<< " mol" <<std::endl;
  //}
  
  return fail;

}


/* *******************************************************************
* Interpolate linearly in time between two values v0 and v1. The time 
* is measuared relative to value v0; so that v1 is at time dt. The
* interpolated data are at time dt_int.            
******************************************************************* */
void CoupledTransport_PK::InterpolateCellVector(
    const Epetra_MultiVector& v0, const Epetra_MultiVector& v1, 
    double dt_int, double dt, Epetra_MultiVector& v_int) 
{
  double a = dt_int / dt;
  double b = 1.0 - a;
  v_int.Update(b, v0, a, v1, 0.);
}

}  // namespace Amanzi
