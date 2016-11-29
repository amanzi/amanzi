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
    using_sorption_isotherms_(false) {};

/* ******************************************************************
* Register fields and evaluators with the State
******************************************************************* */
void Chemistry_PK::Setup(const Teuchos::Ptr<State>& S)
{
  // Require data from flow
  if (!S->HasField("porosity")) {
    S->RequireField("porosity", passwd_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("saturation_liquid")) {
    S->RequireField("saturation_liquid", passwd_)->SetMesh(mesh_)->SetGhosted(false)
      //->SetComponent("cell", AmanziMesh::CELL, 1);
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  }
  
  if (!S->HasField("fluid_density")) {
    S->RequireScalar("fluid_density", passwd_);
  }

  // require transport fields
  std::vector<std::string>::const_iterator it;
  if (!S->HasField("total_component_concentration")) {    
    // set the names for vis
    std::vector<std::vector<std::string> > conc_names_cv(1);
    for (it = comp_names_.begin(); it != comp_names_.end(); ++it) {
      conc_names_cv[0].push_back(*it + std::string(" conc"));
    }
    S->RequireField("total_component_concentration", passwd_, conc_names_cv)
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
    S->RequireField("mineral_volume_fractions", passwd_, vf_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_minerals_);

    S->RequireField("mineral_specific_surface_area", passwd_, ssa_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_minerals_);

    S->RequireField("mineral_rate_constant", passwd_, mrc_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_minerals_);
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
    S->RequireField("sorption_sites", passwd_, ss_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_sorption_sites_);
    S->RequireField("surface_complex_free_site_conc", passwd_, scfsc_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_sorption_sites_);
  }

  if (using_sorption_) {
    S->RequireField("total_sorbed", passwd_)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);

    if (using_sorption_isotherms_) {
      S->RequireField("isotherm_kd", passwd_)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);

      S->RequireField("isotherm_freundlich_n", passwd_)
        ->SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);

      S->RequireField("isotherm_langmuir_b", passwd_)
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

    S->RequireField("free_ion_species", passwd_, species_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);

    S->RequireField("primary_activity_coeff", passwd_, species_names_cv)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_aqueous_components_);
  }

  if (number_ion_exchange_sites_ > 0) {
    S->RequireField("ion_exchange_sites", passwd_)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, number_ion_exchange_sites_);

    S->RequireField("ion_exchange_ref_cation_conc", passwd_)
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
    if (!S->GetField("total_component_concentration", passwd_)->initialized()) {
      InitializeField_("total_component_concentration", 0.0);
    }
    InitializeField_("primary_activity_coeff", 1.0);

    // special initialization of free ion concentration
    if (S_->HasField("free_ion_species")) {
      if (!S_->GetField("free_ion_species")->initialized()) {
        CompositeVector& ion = *S_->GetFieldData("free_ion_species", passwd_);
        const CompositeVector& tcc = *S_->GetFieldData("total_component_concentration");

        ion.Update(0.1, tcc, 0.0);
        S_->GetField("free_ion_species", passwd_)->set_initialized();

        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
          Teuchos::OSTab tab = vo_->getOSTab();
          *vo_->os() << "initilized \"free_ion_species\" to 10% of \"total_component_concentration\"\n";  
        }
      }
    }
    // InitializeField_("free_ion_species", 0.0);

    // Sorption sites: all will have a site density, but we can default to zero
    if (using_sorption_) {
      InitializeField_("total_sorbed", 0.0);
    }

    // Sorption isotherms: Kd required, Langmuir and Freundlich optional
    if (using_sorption_isotherms_) {
      InitializeField_("isotherm_kd", -1.0);
      InitializeField_("isotherm_freundlich_n", 1.0);
      InitializeField_("isotherm_langmuir_b", 1.0);
    }
  }

  // Minerals: vol frac and surface areas
  if (number_minerals_ > 0) {
    InitializeField_("mineral_volume_fractions", 0.0);
    InitializeField_("mineral_specific_surface_area", 1.0);
  }

  // Ion exchange sites: default to 1
  if (number_ion_exchange_sites_ > 0) {
    InitializeField_("ion_exchange_sites", 1.0);
    InitializeField_("ion_exchange_ref_cation_conc", 1.0);
  }

  if (number_sorption_sites_ > 0) {
    InitializeField_("sorption_sites", 1.0);
    InitializeField_("surface_complex_free_site_conc", 1.0);
  }

  // auxiliary fields
  InitializeField_("alquimia_aux_data", 0.0);

  // miscaleneous controls
  initial_conditions_time_ = cp_list_->get<double>("initial conditions time", S->time());
}


/* ******************************************************************
* Process names of materials 
******************************************************************* */
void Chemistry_PK::InitializeField_(std::string fieldname, double default_val)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  if (S_->HasField(fieldname)) {
    if (!S_->GetField(fieldname)->initialized()) {
      S_->GetFieldData(fieldname, passwd_)->PutScalar(default_val);
      S_->GetField(fieldname, passwd_)->set_initialized();
      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
         *vo_->os() << "initilized \"" << fieldname << "\" to value " << default_val << std::endl;  
    }
  }
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
                                           Teuchos::RCP<Teuchos::ParameterList> state_list)
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

  if (state_list->sublist("initial conditions").isSublist("ion_exchange_sites")) {
    // there is currently only at most one site...
    using_sorption_ = true;
    number_ion_exchange_sites_ = 1;
  }

  if (state_list->sublist("initial conditions").isSublist("isotherm_kd")) {
    using_sorption_ = true;
    using_sorption_isotherms_ = true;
  }

  if (state_list->sublist("initial conditions").isSublist("sorption_sites")) {
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

}  // namespace AmanziChemistry
}  // namespace Amanzi
