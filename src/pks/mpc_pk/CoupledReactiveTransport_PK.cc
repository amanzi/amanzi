/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  Process kernel for coupling of Transport PK and Chemistry PK.
*/

#include "CoupledReactiveTransport_PK.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
CoupledReactiveTransport_PK::CoupledReactiveTransport_PK(Teuchos::ParameterList& pk_tree,
                                          const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                          const Teuchos::RCP<State>& S,
                                          const Teuchos::RCP<TreeVector>& soln) :
  ReactiveTransport_PK(pk_tree, global_list, S, soln) { 

  // storage_created = false;
  // chem_step_succeeded = true;

  tranport_pk_ = Teuchos::rcp_dynamic_cast<PK_MPC>(sub_pks_[1]);
  AMANZI_ASSERT(tranport_pk_ != Teuchos::null);

  chemistry_pk_ = Teuchos::rcp_dynamic_cast<PK_MPC>(sub_pks_[0]);
  AMANZI_ASSERT(chemistry_pk_ != Teuchos::null);


}


void CoupledReactiveTransport_PK::Setup(const Teuchos::Ptr<State>& S){

  Amanzi::PK_MPCAdditive<PK>::Setup(S);
  
  cast_sub_pks_();
 
// communicate chemistry engine to transport.
#ifdef ALQUIMIA_ENABLED
  for (int i=0; i<num_sub_domains_; i++){
    tranport_pk_sub_[i]->SetupAlquimia(Teuchos::rcp_static_cast<AmanziChemistry::Alquimia_PK>(chemistry_pk_sub_[i]),
                                       chemistry_pk_sub_[i]->chem_engine());
  }
#endif

}

void CoupledReactiveTransport_PK::cast_sub_pks_(){


  AMANZI_ASSERT(tranport_pk_->num_sub_pks() == chemistry_pk_->num_sub_pks());
  num_sub_domains_ = tranport_pk_->num_sub_pks();

  tranport_pk_sub_.resize(num_sub_domains_);
  chemistry_pk_sub_.resize(num_sub_domains_);

  for (int i=0; i<num_sub_domains_; i++){
    tranport_pk_sub_[i] = Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(tranport_pk_->get_subpk(i));
    AMANZI_ASSERT(tranport_pk_sub_[i] != Teuchos::null);
    chemistry_pk_sub_[i] =  Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(chemistry_pk_->get_subpk(i));
    AMANZI_ASSERT(chemistry_pk_sub_[i] != Teuchos::null);
    //AMANZI_ASSERT(tranport_pk_sub_[i]->domain_name() == chemistry_pk_sub_[i]->domain_name());
  }

  
    

}

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void CoupledReactiveTransport_PK::Initialize(const Teuchos::Ptr<State>& S) {

  Amanzi::PK_MPCAdditive<PK>::Initialize(S);

  for (int i=0; i<num_sub_domains_; i++){
    chemistry_pk_sub_[i] -> Initialize(S);
    tranport_pk_sub_[i] -> Initialize(S);
  }
  
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double CoupledReactiveTransport_PK::get_dt() {

  dTtran_ = tranport_pk_->get_dt();
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
void CoupledReactiveTransport_PK::set_dt(double dt) {

  dTtran_ = dt;
  dTchem_ = dt;
  //dTchem_ = chemistry_pk_->get_dt();
  //if (dTchem_ > dTtran_) dTchem_ = dTtran_;
  //if (dTtran_ > dTchem_) dTtran_= dTchem_; 
  chemistry_pk_->set_dt(dTchem_);
  tranport_pk_->set_dt(dTtran_);
  
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool CoupledReactiveTransport_PK::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;

    
  return fail;
};


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void CoupledReactiveTransport_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {
  chemistry_pk_->CommitStep(t_old, t_new, S);
}



bool CoupledReactiveTransport_PK::AdvanceChemistry_(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
                                                const Epetra_MultiVector& mol_dens,
                                                Teuchos::RCP<Epetra_MultiVector> tcc_copy,
                                                double t_old, double t_new, bool reinit){

    bool pk_fail = false;
  
    chem_pk->set_aqueous_components(tcc_copy);
    //Teuchos::RCP<Teuchos::TimeMonitor> local_flow_monitor = Teuchos::rcp(new Teuchos::TimeMonitor(*alquimia_timer_));
    pk_fail = chem_pk->AdvanceStep(t_old, t_new, reinit);
    //local_flow_monitor = Teuchos::null;    
    *tcc_copy = *chem_pk->aqueous_components();
    
    return pk_fail;

}
  
}  // namespace Amanzi
