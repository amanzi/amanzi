/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Base class for chemical process kernels.
*/
 
// Amanzi
#include "message.hh"

// Chemistry
#include "Chemistry_PK.hh"

namespace Amanzi {
namespace AmanziChemistry {

/* ******************************************************************
* Default constructor that initializes all pointers to NULL
****************************************************************** */
Chemistry_PK::Chemistry_PK() :
    passwd_("state"),
    number_minerals_(0),
    number_ion_exchange_sites_(0),
    number_sorption_sites_(0),
    using_sorption_(false),
    using_sorption_isotherms_(false),    
    number_aqueous_kinetics_(0)
    {};


/* ******************************************************************
* Register fields and evaluators with the State
******************************************************************* */
void Chemistry_PK::Setup(const Teuchos::Ptr<State>& S)
{
  mesh_ = S->GetMesh(domain_);

  saturation_tolerance_ = cp_list_->get<double>("saturation tolerance", 1e-14);
  
  // Require data from flow
  S->RequireField(poro_key_, passwd_)->SetMesh(mesh_)->SetGhosted(false)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(poro_key_);

  S->RequireField(saturation_key_, passwd_)->SetMesh(mesh_)->SetGhosted(false)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(saturation_key_);

  S->RequireField(fluid_den_key_, passwd_)->SetMesh(mesh_)->SetGhosted(false)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(fluid_den_key_);

  // require transport fields
  std::vector<std::string>::const_iterator it;
  if (!S->HasField(tcc_key_)) {    
    // set the names for vis
    std::vector<std::vector<std::string> > conc_names_cv(1);
    for (it = comp_names_.begin(); it != comp_names_.end(); ++it) {
      conc_names_cv[0].push_back(*it + std::string(" conc"));
    }
    S->RequireField(tcc_key_, passwd_, conc_names_cv)
      ->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);
  }

  // require minerals
  if (number_minerals_ > 0) {
    // -- set the names for vis
    std::vector<std::vector<std::string> > vf_names_cv(1);
    std::vector<std::vector<std::string> > ssa_names_cv(1);
    std::vector<std::vector<std::string> > mrc_names_cv(1);

    for (it = mineral_names_.begin(); it != mineral_names_.end(); ++it) {
      vf_names_cv[0].push_back(*it + std::string(" vol frac"));
      ssa_names_cv[0].push_back(*it + std::string(" spec surf area"));
      mrc_names_cv[0].push_back(*it + std::string(" min rate cnst"));
    }

    // -- register two fields
    S->RequireField(min_vol_frac_key_, passwd_, vf_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_minerals_);

    S->RequireField(min_ssa_key_, passwd_, ssa_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_minerals_);

    S->RequireField(mineral_rate_constant_key_, passwd_, mrc_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_minerals_);
  }
  
  // require aqueous kinetics
  if (number_aqueous_kinetics_ > 0) {
    // -- set the names for vis
    std::vector<std::vector<std::string> > aqueous_kinetics_names_cv(1);

    for (it = aqueous_kinetics_names_.begin(); it != aqueous_kinetics_names_.end(); ++it) {
      aqueous_kinetics_names_cv[0].push_back(*it + std::string(" aq kin rate cnst"));
    }

    // -- register field
    S->RequireField(first_order_decay_constant_key_, passwd_, aqueous_kinetics_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_kinetics_);
  }
  
  // require sorption sites
  if (number_sorption_sites_ > 0) {
    // -- set the names for vis
    std::vector<std::vector<std::string> > ss_names_cv(1);
    std::vector<std::vector<std::string> > scfsc_names_cv(1);
    for (it = sorption_site_names_.begin(); it != sorption_site_names_.end(); ++it) {
      ss_names_cv[0].push_back(*it + std::string(" sorption site"));
      scfsc_names_cv[0].push_back(*it + std::string(" surface complex free site conc"));
    }

    // -- register two fields
    S->RequireField(sorp_sites_key_, passwd_, ss_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_sorption_sites_);
    S->RequireField(surf_cfsc_key_, passwd_, scfsc_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_sorption_sites_);
  }

  if (using_sorption_) {
    S->RequireField(total_sorbed_key_, passwd_)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);

    if (using_sorption_isotherms_) {
      S->RequireField(isotherm_kd_key_, passwd_)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);

      S->RequireField(isotherm_freundlich_n_key_, passwd_)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);

      S->RequireField(isotherm_langmuir_b_key_, passwd_)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);
    }
  }

  // ion
  if (number_aqueous_components_ > 0) {
    std::vector<std::vector<std::string> > species_names_cv(1);
    for (it = comp_names_.begin(); it != comp_names_.end(); ++it) {
      species_names_cv[0].push_back(*it);
    }

    S->RequireField(free_ion_species_key_, passwd_, species_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);

    S->RequireField(primary_activity_coeff_key_, passwd_, species_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);
  }

  if (number_ion_exchange_sites_ > 0) {
    S->RequireField(ion_exchange_sites_key_, passwd_)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_ion_exchange_sites_);

    S->RequireField(ion_exchange_ref_cation_conc_key_, passwd_)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_ion_exchange_sites_);
  }
}


/* ******************************************************************
* Most things are initialized through State, but State can only manage that
* if they are always initialized.  If sane defaults are available, or they
* can be derived from other initialized quantities, they are initialized
* here, where we can manage that logic.
******************************************************************* */
void Chemistry_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  // Aqueous species 
  if (number_aqueous_components_ > 0) {
    if (!S->GetField(tcc_key_, passwd_)->initialized()) {
      InitializeField_(S, passwd_, tcc_key_, 0.0);
    }
 
    InitializeField_(S, passwd_, primary_activity_coeff_key_, 1.0);    

    // special initialization of free ion concentration
    if (S_->HasField(free_ion_species_key_)) {
      if (!S_->GetField(free_ion_species_key_)->initialized()) {
        CompositeVector& ion = *S_->GetFieldData(free_ion_species_key_, passwd_);
        const CompositeVector& tcc = *S_->GetFieldData(tcc_key_);

        ion.Update(0.1, tcc, 0.0);
        S_->GetField(free_ion_species_key_, passwd_)->set_initialized();

        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
          Teuchos::OSTab tab = vo_->getOSTab();
          *vo_->os() << "initialized \" "<<free_ion_species_key_<<"\" to 10% of \" "<<
            tcc_key_<<"\"\n";  
        }
      }
    }

    // Sorption sites: all will have a site density, but we can default to zero
    if (using_sorption_) {
      // InitializeField_(S,  total_sorbed_key_, 0.0);
      InitializeField_(S, passwd_, total_sorbed_key_, 0.0);
    }

    // Sorption isotherms: Kd required, Langmuir and Freundlich optional
    if (using_sorption_isotherms_) {
      InitializeField_(S, passwd_, isotherm_kd_key_, -1.0);
      InitializeField_(S, passwd_, isotherm_freundlich_n_key_, 1.0);
      InitializeField_(S, passwd_, isotherm_langmuir_b_key_, 1.0);
    }
  }

  // Minerals: vol frac and surface areas
  if (number_minerals_ > 0) {
    InitializeField_(S, passwd_, min_vol_frac_key_, 0.0);
    InitializeField_(S, passwd_, min_ssa_key_, 1.0);

  }

  // Aqueous kinetics
  if (number_aqueous_kinetics_ > 0) {
    InitializeField_(S, passwd_, first_order_decay_constant_key_, 0.0);
  }
  
  // Ion exchange sites: default to 1
  if (number_ion_exchange_sites_ > 0) {
    InitializeField_(S, passwd_, ion_exchange_sites_key_, 1.0);
    InitializeField_(S, passwd_, ion_exchange_ref_cation_conc_key_, 1.0);    
  }

  if (number_sorption_sites_ > 0) {
    InitializeField_(S, passwd_, sorp_sites_key_, 1.0);
    InitializeField_(S, passwd_, surf_cfsc_key_, 1.0);
  }

  // auxiliary fields
  InitializeField_(S, passwd_, alquimia_aux_data_key_, 0.0);

  // miscaleneous controls
  initial_conditions_time_ = cp_list_->get<double>("initial conditions time", S->time());
}


/* ******************************************************************
* Process names of materials 
******************************************************************* */
void Chemistry_PK::InitializeMinerals(Teuchos::RCP<Teuchos::ParameterList> plist)
{
  mineral_names_.clear();
  if (plist->isParameter("minerals")) {
    mineral_names_ = plist->get<Teuchos::Array<std::string> >("minerals").toVector();
  }

  number_minerals_ = mineral_names_.size();
}


/* ******************************************************************
* Process names of sorption sites
* NOTE: Do we need to worry about sorption sites?
******************************************************************* */
void Chemistry_PK::InitializeSorptionSites(Teuchos::RCP<Teuchos::ParameterList> plist,
                                           Teuchos::RCP<Teuchos::ParameterList> ic_list)
{
  sorption_site_names_.clear();
  if (plist->isParameter("sorption sites")) {
    sorption_site_names_ = plist->get<Teuchos::Array<std::string> >("sorption sites").toVector();
  }

  number_sorption_sites_ = sorption_site_names_.size();
  using_sorption_ = (number_sorption_sites_ > 0);

  // check if there is an initial condition for ion_exchange_sites
  number_ion_exchange_sites_ = 0;
  using_sorption_isotherms_ = false;

  if (ic_list->isSublist(ion_exchange_sites_key_)) {
    // there is currently only at most one site...
    using_sorption_ = true;
    number_ion_exchange_sites_ = 1;
  }

  if (ic_list->isSublist(isotherm_kd_key_)) {
    using_sorption_ = true;
    using_sorption_isotherms_ = true;
  }

  if (ic_list->isSublist(sorp_sites_key_)) {
    using_sorption_ = true;
  }

  // in the old version, this was only in the Block sublist... may need work?
  if (plist->isParameter("Cation Exchange Capacity")) {
    using_sorption_ = true;
    number_ion_exchange_sites_ = 1;
  }
}


/* *******************************************************************
* I/O or error messages
******************************************************************* */
void Chemistry_PK::ErrorAnalysis(int ierr, std::string& internal_msg)
{
  int tmp_out[2], tmp_in[2] = {ierr, (int) internal_msg.size()};
  mesh_->get_comm()->MaxAll(tmp_in, tmp_out, 2);

  if (tmp_out[0] != 0) {
    // update time control parameters
    num_successful_steps_ = 0;

    // get at least one error message
    int n = tmp_out[1];
    int msg_out[n + 1], msg_in[n + 1], m(mesh_->get_comm()->MyPID());
    internal_msg.resize(n);

    Errors::encode_string(internal_msg, n, m, msg_in);
    mesh_->get_comm()->MaxAll(msg_in, msg_out, n + 1);
    Errors::decode_string(msg_out, n, internal_msg);

    vo_->Write(Teuchos::VERB_HIGH, internal_msg);

    Errors::Message msg(internal_msg);
    Exceptions::amanzi_throw(msg); 
  }  
}

void Chemistry_PK::CopyFieldstoNewState(const Teuchos::RCP<State>& S_next) {

  if (S_->HasField(min_vol_frac_key_)&& S_next->HasField(min_vol_frac_key_)) {
    *S_next->GetFieldData(min_vol_frac_key_, passwd_)->ViewComponent("cell", false) =
        *S_->GetFieldData(min_vol_frac_key_, passwd_)->ViewComponent("cell", false);
  }

  if (S_->HasField(min_ssa_key_)&& S_next->HasField(min_ssa_key_)) {
    *S_next->GetFieldData(min_ssa_key_, passwd_)->ViewComponent("cell", false) =
        *S_->GetFieldData(min_ssa_key_, passwd_)->ViewComponent("cell", false);
  }

  if (S_->HasField(sorp_sites_key_)&& S_next->HasField(sorp_sites_key_)) {
    *S_next->GetFieldData(sorp_sites_key_, passwd_)->ViewComponent("cell", false) =
        *S_->GetFieldData(sorp_sites_key_, passwd_)->ViewComponent("cell", false);
  }

  if (S_->HasField(surf_cfsc_key_)&& S_next->HasField(surf_cfsc_key_)) {
    *S_next->GetFieldData(surf_cfsc_key_, passwd_)->ViewComponent("cell", false) =
      *S_->GetFieldData(surf_cfsc_key_, passwd_)->ViewComponent("cell", false);
  }

  if (S_->HasField(total_sorbed_key_)&& S_next->HasField(total_sorbed_key_)) {
    *S_next->GetFieldData(total_sorbed_key_, passwd_)->ViewComponent("cell", false) =
      *S_->GetFieldData(total_sorbed_key_, passwd_)->ViewComponent("cell", false);
  }

  if (S_->HasField(isotherm_kd_key_)&& S_next->HasField(isotherm_kd_key_)) {
    *S_next->GetFieldData(isotherm_kd_key_, passwd_)->ViewComponent("cell", false) =
      *S_->GetFieldData(isotherm_kd_key_, passwd_)->ViewComponent("cell", false);
  }

  if (S_->HasField(isotherm_freundlich_n_key_)&& S_next->HasField(isotherm_freundlich_n_key_)) {
    *S_next->GetFieldData(isotherm_freundlich_n_key_, passwd_)->ViewComponent("cell", false) =
      *S_->GetFieldData(isotherm_freundlich_n_key_, passwd_)->ViewComponent("cell", false);
  }
    
  if (S_->HasField(isotherm_langmuir_b_key_)&& S_next->HasField(isotherm_langmuir_b_key_)) {
    *S_next->GetFieldData(isotherm_langmuir_b_key_, passwd_)->ViewComponent("cell", false) =
      *S_->GetFieldData(isotherm_langmuir_b_key_, passwd_)->ViewComponent("cell", false);
  }

  if (S_->HasField(free_ion_species_key_)&& S_next->HasField(free_ion_species_key_)) {
    *S_next->GetFieldData(free_ion_species_key_, passwd_)->ViewComponent("cell", false) =
      *S_->GetFieldData(free_ion_species_key_, passwd_)->ViewComponent("cell", false);
  }

  if (S_->HasField(primary_activity_coeff_key_)&& S_next->HasField(primary_activity_coeff_key_)) {
    *S_next->GetFieldData(primary_activity_coeff_key_, passwd_)->ViewComponent("cell", false) =
      *S_->GetFieldData(primary_activity_coeff_key_, passwd_)->ViewComponent("cell", false);
  }
    
  if (S_->HasField(ion_exchange_sites_key_)&& S_next->HasField(ion_exchange_sites_key_)) {
    *S_next->GetFieldData(ion_exchange_sites_key_, passwd_)->ViewComponent("cell", false) =
      *S_->GetFieldData(ion_exchange_sites_key_, passwd_)->ViewComponent("cell", false);
  }

  if (S_->HasField(ion_exchange_ref_cation_conc_key_)&& S_next->HasField(ion_exchange_ref_cation_conc_key_)) {
    *S_next->GetFieldData(ion_exchange_ref_cation_conc_key_, passwd_)->ViewComponent("cell", false) =
      *S_->GetFieldData(ion_exchange_ref_cation_conc_key_, passwd_)->ViewComponent("cell", false);
  }


 

}

}  // namespace AmanziChemistry
}  // namespace Amanzi
