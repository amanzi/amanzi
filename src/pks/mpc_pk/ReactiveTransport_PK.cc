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
  Amanzi::PK_MPCAdditive<PK>(pk_tree, global_list, S, soln) { 

  storage_created = false;
  chem_step_succeeded = true;

  transport_pk_ = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(sub_pks_[1]);
  AMANZI_ASSERT(transport_pk_ != Teuchos::null);

  chemistry_pk_ = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(sub_pks_[0]);
  AMANZI_ASSERT(chemistry_pk_ != Teuchos::null);

  // communicate chemistry engine to transport.
#ifdef ALQUIMIA_ENABLED
  transport_pk_->SetupAlquimia(Teuchos::rcp_static_cast<AmanziChemistry::Alquimia_PK>(chemistry_pk_),
                              chemistry_pk_->chem_engine());
#endif

  // master_ = 1;  // Transport;
  // slave_ = 0;  // Chemistry;
}


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void ReactiveTransport_PK::Initialize(const Teuchos::Ptr<State>& S) {
  Amanzi::PK_MPCAdditive<PK>::Initialize(S);

  if (S->HasField("total_component_concentration")) {
    total_component_concentration_stor = 
       Teuchos::rcp(new Epetra_MultiVector(*S->GetFieldData("total_component_concentration")
                                             ->ViewComponent("cell", true)));
    storage_created = true;
  }
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double ReactiveTransport_PK::get_dt() {

  dTtran_ = transport_pk_->get_dt();
  dTchem_ = chemistry_pk_->get_dt();

  if (!chem_step_succeeded && (dTchem_/dTtran_ > 0.99)) {
     dTchem_ *= 0.5;
  } 

  if (dTtran_ > dTchem_) dTtran_= dTchem_; 
  
  return dTchem_;
}


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void ReactiveTransport_PK::set_dt(double dt) {
  dTtran_ = dt;
  dTchem_ = dt;
  //dTchem_ = chemistry_pk_->get_dt();
  //if (dTchem_ > dTtran_) dTchem_ = dTtran_;
  //if (dTtran_ > dTchem_) dTtran_= dTchem_; 
  chemistry_pk_->set_dt(dTchem_);
  transport_pk_->set_dt(dTtran_);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool ReactiveTransport_PK::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;
  chem_step_succeeded = false;

  // First we do a transport step.
  bool pk_fail = transport_pk_->AdvanceStep(t_old, t_new, reinit);

  // Right now transport step is always succeeded.
  if (!pk_fail) {
    *total_component_concentration_stor = *transport_pk_->total_component_concentration()->ViewComponent("cell", true);
  } else {
    Errors::Message message("MPC: Transport PK returned an unexpected error.");
    Exceptions::amanzi_throw(message);
  }

  // Second, we do a chemistry step.
  try {
    chemistry_pk_->set_aqueous_components(total_component_concentration_stor);

    pk_fail = chemistry_pk_->AdvanceStep(t_old, t_new, reinit);
    chem_step_succeeded = true;
 
    *S_->GetFieldData("total_component_concentration", "state")
       ->ViewComponent("cell", true) = *chemistry_pk_->aqueous_components();
  }
  catch (const Errors::Message& chem_error) {
    fail = true;
  }
    
  return fail;
};


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void ReactiveTransport_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {
  chemistry_pk_->CommitStep(t_old, t_new, S);
}

}  // namespace Amanzi

