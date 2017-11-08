/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  Process kernel for coupling of Transport PK and Chemistry PK.
*/

#include "mpc_reactivetransport_pk.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
ReactiveTransport_PK_ATS::ReactiveTransport_PK_ATS(Teuchos::ParameterList& pk_tree,
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
  rt_pk_list_ = Teuchos::sublist(pk_list, pk_name, true);
  
  transport_pk_index_ = rt_pk_list_->get<int>("transport index", 0);
  chemistry_pk_index_ = rt_pk_list_->get<int>("chemistry index", 1 - transport_pk_index_);

  tranport_pk_ = Teuchos::rcp_dynamic_cast<Transport::Transport_PK_ATS>(sub_pks_[transport_pk_index_]);
  ASSERT(tranport_pk_ != Teuchos::null);

  chemistry_pk_ = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(sub_pks_[chemistry_pk_index_]);
  ASSERT(chemistry_pk_ != Teuchos::null);

  // communicate chemistry engine to transport.
#ifdef ALQUIMIA_ENABLED
  tranport_pk_->SetupAlquimia(Teuchos::rcp_static_cast<AmanziChemistry::Alquimia_PK>(chemistry_pk_),
                              chemistry_pk_->chem_engine());
#endif

  // std::cout<<tranport_pk_->get_domain_name()<<"\n";
  // std::cout<<chemistry_pk_->get_domain_name()<<"\n";

  // ASSERT(tranport_pk_->get_domain_name() == chemistry_pk_->get_domain_name());
  // master_ = 1;  // Transport;
  // slave_ = 0;  // Chemistry;
}


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void ReactiveTransport_PK_ATS::Initialize(const Teuchos::Ptr<State>& S) {
  Amanzi::PK_MPCAdditive<PK>::Initialize(S);

  // if (S->HasField("total_component_concentration")) {
  //   total_component_concentration_stor = 
  //      Teuchos::rcp(new Epetra_MultiVector(*S->GetFieldData("total_component_concentration")
  //                                            ->ViewComponent("cell", true)));
  //   storage_created = true;
  // }
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double ReactiveTransport_PK_ATS::get_dt() {

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
void ReactiveTransport_PK_ATS::set_dt(double dt) {
  dTtran_ = dt;
  dTchem_ = dt;
  //dTchem_ = chemistry_pk_ -> get_dt();
  //if (dTchem_ > dTtran_) dTchem_ = dTtran_;
  //if (dTtran_ > dTchem_) dTtran_= dTchem_; 
  chemistry_pk_->set_dt(dTchem_);
  tranport_pk_->set_dt(dTtran_);
}

void ReactiveTransport_PK_ATS::set_states(const Teuchos::RCP<const State>& S,
                                          const Teuchos::RCP<State>& S_inter,
                                          const Teuchos::RCP<State>& S_next) {
  //  PKDefaultBase::set_states(S, S_inter, S_next);
  //S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;

  //chemistry_pk_->set_states(S, S_inter, S_next);
  tranport_pk_->set_states(S, S_inter, S_next);

}

// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool ReactiveTransport_PK_ATS::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;
  chem_step_succeeded = false;

  Key domain_name = tranport_pk_->get_domain_name();
  Key tcc_key = Keys::getKey(domain_name, "total_component_concentration");
  Key mol_den_key = Keys::getKey(domain_name,  "molar_density_liquid");

  // First we do a transport step.
  bool pk_fail = tranport_pk_->AdvanceStep(t_old, t_new, reinit);

  if (pk_fail){
    Errors::Message message("MPC: Transport PK returned an unexpected error.");
    Exceptions::amanzi_throw(message);
  } 

  // Second, we do a chemistry step.
  try {
    Teuchos::RCP<Epetra_MultiVector> tcc_copy =
      S_->GetFieldCopyData(tcc_key,"subcycling","state")->ViewComponent("cell", true);

    Teuchos::RCP<const Epetra_MultiVector> mol_dens =
      S_->GetFieldData(mol_den_key)->ViewComponent("cell", true);

    Teuchos::RCP<const AmanziMesh::Mesh> mesh = S_->GetMesh(domain_name);
    int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

    // convert from mole fraction[-] to mol/L
    for (int c=0; c<ncells_owned; c++) (*tcc_copy)[0][c] *= (*mol_dens)[0][c] / 1000.;

    
    chemistry_pk_->set_aqueous_components(tcc_copy);
    pk_fail = chemistry_pk_->AdvanceStep(t_old, t_new, reinit);
    chem_step_succeeded = true;  
 
    *tcc_copy = *chemistry_pk_->aqueous_components();
    // convert from mol/L fraction to mole fraction[-]
    for (int c=0; c<ncells_owned; c++) (*tcc_copy)[0][c] /= (*mol_dens)[0][c] / 1000.;
    
  }
  catch (const Errors::Message& chem_error) {
    fail = true;
  }
    
  return fail;
};


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void ReactiveTransport_PK_ATS::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {

  tranport_pk_->CommitStep(t_old, t_new, S);
  chemistry_pk_->CommitStep(t_old, t_new, S);
}

}  // namespace Amanzi

