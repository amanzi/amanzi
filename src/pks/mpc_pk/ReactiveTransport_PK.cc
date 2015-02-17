/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  Process kernel for coupling of Transport PK and Chemistry PK.
*/

#include "ReactiveTransport_PK.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
ReactiveTransport_PK::ReactiveTransport_PK(Teuchos::ParameterList& pk_tree,
					   const Teuchos::RCP<Teuchos::ParameterList>& global_list,
					   const Teuchos::RCP<State>& S,
					   const Teuchos::RCP<TreeVector>& soln) :
    Amanzi::MPCWeak(pk_tree, global_list, S, soln) { 

  storage_created = false;
  chem_step_succeeded = true;

  tranport_pk_ = Teuchos::rcp_dynamic_cast<Transport::Transport_PK_Wrapper>(sub_pks_[1]);
  ASSERT(tranport_pk_ != Teuchos::null);
  chemistry_pk_ = Teuchos::rcp_dynamic_cast<Amanzi::AmanziChemistry::Chemistry_PK_Wrapper>(sub_pks_[0]);
  ASSERT(chemistry_pk_ != Teuchos::null);

  // master_ = 1; // Transport;
  // slave_ = 0; // Chemistry;
}


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void ReactiveTransport_PK::Initialize() {

  Amanzi::MPCWeak::Initialize();

  if (S_->HasField("total_component_concentration")) {
    total_component_concentration_stor = 
       Teuchos::rcp(new Epetra_MultiVector(*S_->GetFieldData("total_component_concentration")->ViewComponent("cell", true)));
      storage_created  = true;
  }
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double ReactiveTransport_PK::get_dt() {

  dTtran_ = tranport_pk_ -> get_dt();
  dTchem_ = chemistry_pk_ -> get_dt();

  if (!chem_step_succeeded && (dTchem_/dTtran_ > 0.99)) {
     dTchem_ *= 0.5;
  } 

  //std::cout<<"Transport dT="<<dTtran_<<" Chemistry dT="<<dTchem_<<"\n";
  if (dTtran_ > dTchem_) dTtran_= dTchem_; 
  
  return dTchem_;
};


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void ReactiveTransport_PK::set_dt(double dt) {
  dTtran_ = dt;
  dTchem_ = dt;
  //dTchem_ = chemistry_pk_ -> get_dt();
  //if (dTchem_ > dTtran_) dTchem_ = dTtran_;
  //if (dTtran_ > dTchem_) dTtran_= dTchem_; 
  chemistry_pk_ -> set_dt(dTchem_);
  tranport_pk_ -> set_dt(dTtran_);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool ReactiveTransport_PK::AdvanceStep(double t_old, double t_new) {
  bool fail = false;
  chem_step_succeeded = false;

  int ok = tranport_pk_->AdvanceStep(t_old, t_new);

  if (ok == 0) {
    *total_component_concentration_stor = *tranport_pk_->total_component_concentration()->ViewComponent("cell", true);
  } else {
    Errors::Message message("MPC: error... Transport_PK.advance returned an error status");
    Exceptions::amanzi_throw(message);
  }
  chemistry_pk_->set_total_component_concentration(total_component_concentration_stor);

  ok = chemistry_pk_->AdvanceStep(t_old, t_new);
  chem_step_succeeded = true;
    
  return fail;
};

}  // namespace Amanzi

