/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

*/

#include "mpc_morphology_pk.hh"


namespace Amanzi {

  Morphology_PK::Morphology_PK(Teuchos::ParameterList& pk_tree_or_fe_list,
                               const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                               const Teuchos::RCP<State>& S,
                               const Teuchos::RCP<TreeVector>& soln) :
    PK(pk_tree_or_fe_list, global_list, S, soln),
    WeakMPC(pk_tree_or_fe_list, global_list, S, soln)
  {

     // Create verbosity object.
    vo_ = Teuchos::null;
    Teuchos::ParameterList vlist;
    vlist.sublist("verbose object") = global_list -> sublist("verbose object");
    vo_ =  Teuchos::rcp(new VerboseObject("Morphology_PK", vlist)); 
    
    name_ = "morphology pk";

    Teuchos::Array<std::string> pk_order = plist_->get< Teuchos::Array<std::string> >("PKs order");
  
  }

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double Morphology_PK::get_dt() {

  
  double flow_dt = sub_pks_[0]->get_dt() ;
  set_dt(flow_dt);
  return flow_dt;
 
}


void Morphology_PK::Setup(const Teuchos::Ptr<State>& S){

  //passwd_ = "coupled_transport";  // owner's password
  //passwd_ = "state";  // owner's password

  WeakMPC::Setup(S);
}

void Morphology_PK::Initialize(const Teuchos::Ptr<State>& S){

  WeakMPC::Initialize(S);

}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool Morphology_PK::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;

  double dt_MPC = S_->final_time() - S_->initial_time();

  sub_pks_[0]->AdvanceStep(t_old, t_new, reinit);
  
  return fail;

}


}  // namespace Amanzi
