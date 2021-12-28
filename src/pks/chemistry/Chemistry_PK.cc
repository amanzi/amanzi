/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Base class for chemical process kernels.
*/

#include "message.hh"

// Chemistry
#include "Chemistry_PK.hh"

namespace Amanzi {
namespace AmanziChemistry {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

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
  
  // require flow fields
  S->Require<CV_t, CVS_t>(poro_key_, Tags::DEFAULT, passwd_).SetMesh(mesh_)
    ->SetGhosted(false)->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireEvaluator(poro_key_);

  S->Require<CV_t, CVS_t>(saturation_key_, Tags::DEFAULT, passwd_).SetMesh(mesh_)
    ->SetGhosted(false)->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireEvaluator(saturation_key_);

  S->Require<CV_t, CVS_t>(fluid_den_key_, Tags::DEFAULT, passwd_).SetMesh(mesh_)
    ->SetGhosted(false)->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireEvaluator(fluid_den_key_);

  // require transport fields
  std::vector<std::string>::const_iterator it;
  if (!S->HasData(tcc_key_)) {
    S->Require<CV_t, CVS_t>(tcc_key_, Tags::DEFAULT, passwd_, comp_names_)
      .SetMesh(mesh_)->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);
  }

  // require minerals
  if (number_minerals_ > 0) {
    S->Require<CV_t, CVS_t>(min_vol_frac_key_, Tags::DEFAULT, passwd_, mineral_names_)
      .SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_minerals_);

    S->Require<CV_t, CVS_t>(min_ssa_key_, Tags::DEFAULT, passwd_, mineral_names_)
      .SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_minerals_);

    S->Require<CV_t, CVS_t>(mineral_rate_constant_key_, Tags::DEFAULT, passwd_, mineral_names_)
      .SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_minerals_);
  }
  
  // require aqueous kinetics
  if (number_aqueous_kinetics_ > 0) {
    std::vector<std::string> subfield_names;

    for (it = aqueous_kinetics_names_.begin(); it != aqueous_kinetics_names_.end(); ++it) {
      subfield_names.push_back(*it + std::string(" aq kin rate cnst"));
    }

    // -- register field
    S->Require<CV_t, CVS_t>(first_order_decay_constant_key_, Tags::DEFAULT, passwd_, subfield_names)
      .SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_kinetics_);
  }
  
  // require sorption sites
  if (number_sorption_sites_ > 0) {
    std::vector<std::string> ss_names_cv, scfsc_names_cv;
    for (it = sorption_site_names_.begin(); it != sorption_site_names_.end(); ++it) {
      ss_names_cv.push_back(*it + std::string(" sorption site"));
      scfsc_names_cv.push_back(*it + std::string(" surface complex free site conc"));
    }

    // -- register two fields
    S->Require<CV_t, CVS_t>(sorp_sites_key_, Tags::DEFAULT, passwd_, ss_names_cv)
      .SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_sorption_sites_);
    S->Require<CV_t, CVS_t>(surf_cfsc_key_, Tags::DEFAULT, passwd_, scfsc_names_cv)
      .SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_sorption_sites_);
  }

  if (using_sorption_) {
    S->Require<CV_t, CVS_t>(total_sorbed_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);

    if (using_sorption_isotherms_) {
      S->Require<CV_t, CVS_t>(isotherm_kd_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);

      S->Require<CV_t, CVS_t>(isotherm_freundlich_n_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);

      S->Require<CV_t, CVS_t>(isotherm_langmuir_b_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);
    }
  }

  // ion
  if (number_aqueous_components_ > 0) {
    S->Require<CV_t, CVS_t>(free_ion_species_key_, Tags::DEFAULT, passwd_, comp_names_)
      .SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);

    S->Require<CV_t, CVS_t>(primary_activity_coeff_key_, Tags::DEFAULT, passwd_, comp_names_)
      .SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);
  }

  if (number_ion_exchange_sites_ > 0) {
    S->Require<CV_t, CVS_t>(ion_exchange_sites_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_ion_exchange_sites_);

    S->Require<CV_t, CVS_t>(ion_exchange_ref_cation_conc_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)->SetGhosted(false)
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
    if (!S->GetRecordW(tcc_key_, passwd_).initialized()) {
      InitializeField_(S, passwd_, tcc_key_, 0.0);
    }
 
    InitializeField_(S, passwd_, primary_activity_coeff_key_, 1.0);    
    InitializeField_(S, passwd_, free_ion_species_key_, 1.0e-9);    

    // Sorption sites: all will have a site density, but we can default to zero
    if (using_sorption_) {
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

  KeyTriple split;
  bool is_ds = Keys::splitDomainSet(isotherm_kd_key_, split);
  if (is_ds) {
    Key lifted_key = Keys::getKey(std::get<0>(split), "*", std::get<2>(split));
    
    if (ic_list->isSublist(lifted_key)) {
      using_sorption_ = true;
      using_sorption_isotherms_ = true;
    }
  }
  
  if (ic_list->isSublist(sorp_sites_key_)) {
    using_sorption_ = true;
  }

  // in the old version, this was only in the Block sublist... may need work?
  if (plist->isParameter("Cation Exchange Capacity")) {
    using_sorption_ = true;
    number_ion_exchange_sites_ = 1;
  }

  if (glist_->isSublist("thermodynamic database")) {
    if (glist_->sublist("thermodynamic database").isSublist("isotherms")) {
      using_sorption_ = true;
    }
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


void Chemistry_PK::CopyFieldstoNewState(const Teuchos::RCP<State>& S_next)
{
  std::vector<std::string> fields = {
    min_vol_frac_key_,
    min_ssa_key_,
    sorp_sites_key_,
    surf_cfsc_key_,
    total_sorbed_key_,
    isotherm_kd_key_,
    isotherm_freundlich_n_key_,
    isotherm_langmuir_b_key_,
    free_ion_species_key_,
    primary_activity_coeff_key_,
    ion_exchange_sites_key_,
    ion_exchange_ref_cation_conc_key_
  };

  for (int i = 0; i < fields.size(); ++i) {
    if (S_->HasData(fields[i]) && S_next->HasData(fields[i])) {
      *S_next->GetW<CompositeVector>(fields[i], Tags::DEFAULT, passwd_).ViewComponent("cell") =
        *S_->Get<CompositeVector>(fields[i]).ViewComponent("cell");
    }
  }
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
