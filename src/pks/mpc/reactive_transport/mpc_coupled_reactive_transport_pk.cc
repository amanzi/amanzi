/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  Process kernel for coupling of Transport PK and Chemistry PK.
*/

#include "mpc_coupled_reactive_transport_pk.hh"

namespace Amanzi { 

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
Coupled_ReactiveTransport_PK_ATS::Coupled_ReactiveTransport_PK_ATS(
                                          Teuchos::ParameterList& pk_tree,
                                          const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                          const Teuchos::RCP<State>& S,
                                          const Teuchos::RCP<TreeVector>& soln) :
  Amanzi::PK_MPCAdditive<PK>(pk_tree, global_list, S, soln)
 { 

  storage_created = false;
  chem_step_succeeded = true;
  std::string pk_name = pk_tree.name();
  
  boost::iterator_range<std::string::iterator> res = boost::algorithm::find_last(pk_name,"->"); 
  if (res.end() - pk_name.end() != 0) boost::algorithm::erase_head(pk_name,  res.end() - pk_name.begin());

// Create miscaleneous lists.
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(global_list, "PKs", true);
  //std::cout<<*pk_list;
  crt_pk_list_ = Teuchos::sublist(pk_list, pk_name, true);
  
  transport_pk_index_ = crt_pk_list_->get<int>("transport index", 0);
  chemistry_pk_index_ = crt_pk_list_->get<int>("chemistry index", 1 - transport_pk_index_);

  tranport_pk_ = Teuchos::rcp_dynamic_cast<CoupledTransport_PK>(sub_pks_[transport_pk_index_]);
  ASSERT(tranport_pk_ != Teuchos::null);
  
  chemistry_pk_ = Teuchos::rcp_dynamic_cast<WeakMPC>(sub_pks_[chemistry_pk_index_]);
  ASSERT(chemistry_pk_ != Teuchos::null);

  tranport_pk_subsurface_ = 
    Teuchos::rcp_dynamic_cast<Transport::Transport_PK_ATS>(tranport_pk_->get_subpk(0));
  ASSERT(tranport_pk_subsurface_!= Teuchos::null);
  tranport_pk_overland_ = 
    Teuchos::rcp_dynamic_cast<Transport::Transport_PK_ATS>(tranport_pk_->get_subpk(1));
  ASSERT(tranport_pk_overland_!= Teuchos::null);

  chemistry_pk_subsurface_ = 
    Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(chemistry_pk_->get_subpk(0));
  ASSERT(chemistry_pk_subsurface_!= Teuchos::null);
  chemistry_pk_overland_ = 
    Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(chemistry_pk_->get_subpk(1));
  ASSERT(chemistry_pk_overland_!= Teuchos::null);

  // std::cout<<tranport_pk_subsurface_->name()<<"\n";
  // std::cout<<tranport_pk_overland_->name()<<"\n";
  // std::cout<<chemistry_pk_subsurface_->name()<<"\n";
  // std::cout<<chemistry_pk_overland_->name()<<"\n";


 }

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double Coupled_ReactiveTransport_PK_ATS::get_dt() {

  dTtran_ = tranport_pk_->get_dt();
  dTchem_ = chemistry_pk_->get_dt();

  if (!chem_step_succeeded && (dTchem_/dTtran_ > 0.99)) {
     dTchem_ *= 0.5;
  } 

  if (dTtran_ > dTchem_) dTtran_= dTchem_; 
  
  return dTchem_;
}

void Coupled_ReactiveTransport_PK_ATS::set_states(const Teuchos::RCP<const State>& S,
                                                  const Teuchos::RCP<State>& S_inter,
                                                  const Teuchos::RCP<State>& S_next) {
  //  PKDefaultBase::set_states(S, S_inter, S_next);
  //S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;

  //chemistry_pk_->set_states(S, S_inter, S_next);
  tranport_pk_->set_states(S, S_inter, S_next);

}



void Coupled_ReactiveTransport_PK_ATS::Setup(const Teuchos::Ptr<State>& S){

  Amanzi::PK_MPCAdditive<PK>::Setup(S);

  // communicate chemistry engine to transport.
#ifdef ALQUIMIA_ENABLED
  tranport_pk_subsurface_->SetupAlquimia(Teuchos::rcp_static_cast<AmanziChemistry::Alquimia_PK>(chemistry_pk_subsurface_),
                                            chemistry_pk_subsurface_->chem_engine());
  tranport_pk_overland_->SetupAlquimia(Teuchos::rcp_static_cast<AmanziChemistry::Alquimia_PK>(chemistry_pk_overland_),
                                            chemistry_pk_overland_->chem_engine());
#endif


}


bool Coupled_ReactiveTransport_PK_ATS::AdvanceStep(double t_old, double t_new, bool reinit) {

  bool fail = false;
  chem_step_succeeded = false;
  Key subsurface_domain_key = tranport_pk_subsurface_->get_domain_name();
  Key overland_domain_key = tranport_pk_overland_->get_domain_name();
  Key tcc_sub_key = getKey(subsurface_domain_key, "total_component_concentration");
  Key tcc_over_key = getKey(overland_domain_key, "total_component_concentration");

  // First we do a transport step.
  bool pk_fail = tranport_pk_->AdvanceStep(t_old, t_new, reinit);

  if (pk_fail){
    Errors::Message message("MPC: Coupled Transport PK returned an unexpected error.");
    Exceptions::amanzi_throw(message);
  }


  try {
    Teuchos::RCP<Epetra_MultiVector> tcc_sub = 
      S_->GetFieldCopyData(tcc_sub_key,"subcycling","state")->ViewComponent("cell", true);
    chemistry_pk_subsurface_->set_aqueous_components(tcc_sub);

    Teuchos::RCP<Epetra_MultiVector> tcc_over =
      S_->GetFieldCopyData(tcc_over_key,"subcycling", "state")->ViewComponent("cell", true);
    chemistry_pk_overland_->set_aqueous_components(tcc_over);

    pk_fail = chemistry_pk_->AdvanceStep(t_old, t_new, reinit);
    chem_step_succeeded = true;

    //*S_->GetFieldData(tcc_sub_key, "state")->ViewComponent("cell", true) = *chemistry_pk_subsurface_->aqueous_components();
    *tcc_sub = *chemistry_pk_subsurface_->aqueous_components();

    // *S_->GetFieldData(tcc_over_key, "state")->ViewComponent("cell", true)
    //   = *chemistry_pk_overland_->aqueous_components();
    *tcc_over = *chemistry_pk_overland_->aqueous_components();
    
  }
  catch (const Errors::Message& chem_error) {
        fail = true;
  }
    
  return fail;
};


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void Coupled_ReactiveTransport_PK_ATS::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {

  tranport_pk_->CommitStep(t_old, t_new, S);
  chemistry_pk_->CommitStep(t_old, t_new, S);

}

   

}// namespace
