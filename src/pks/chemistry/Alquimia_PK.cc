/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry PK

  Trilinos based process kernel for chemistry. All geochemistry
  calculations live in the chemistry library. The PK stores the
  instance of the chemistry object and drives the chemistry
  calculations on a cell by cell basis. It handles the movement of
  data back and forth between the amanzi memory and the chemistry
  library data structures.
*/

#include <algorithm>
#include <set>
#include <string>

// TPLs
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "errors.hh"
#include "exceptions.hh"
#include "Mesh.hh"
#include "pk_helpers.hh"

// Chemistry
#include "Alquimia_PK.hh"

namespace Amanzi {
namespace AmanziChemistry {


/* *******************************************************************
* Constructor
******************************************************************* */
Alquimia_PK::Alquimia_PK(Teuchos::ParameterList& pk_tree,
                         const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, global_list, S, soln),
    Chemistry_PK(pk_tree, global_list, S, soln),
    chem_initialized_(false),
    number_aqueous_kinetics_(0),
    number_sorption_sites_(0),
    number_isotherm_species_(0),
    number_ion_exchange_sites_(0)
{
  std::string engine_name = plist_->get<std::string>("engine");
  std::string engine_inputfile = plist_->get<std::string>("engine input file");
  chem_engine_ = Teuchos::rcp(new AmanziChemistry::ChemistryEngine(engine_name, engine_inputfile));
  chem_engine_->InitState(beaker_.properties, beaker_.state, beaker_.aux_data, beaker_.aux_output);

  // Amanzi stores many things in the state->initial conditions list, not locally.  Move them locally.
  if (global_list->sublist("state").sublist("initial conditions").isSublist("geochemical conditions")) {
    plist_->sublist("initial conditions").sublist("geochemical conditions") =
      global_list->sublist("state").sublist("initial conditions").sublist("geochemical conditions");
  }

}


/* *******************************************************************
* Destroy ansilary data structures.
******************************************************************* */
Alquimia_PK::~Alquimia_PK()
{
  chem_engine_->FreeState(beaker_.properties, beaker_.state, beaker_.aux_data, beaker_.aux_output);
}


/* *******************************************************************
* Parser
******************************************************************* */
void
Alquimia_PK::parseParameterList()
{
  Chemistry_PK::parseParameterList();

  // grab the component names and sizes
  chem_engine_->GetPrimarySpeciesNames(aqueous_comp_names_);
  number_aqueous_components_ = aqueous_comp_names_.size();

  chem_engine_->GetMineralNames(mineral_comp_names_);
  number_mineral_components_ = mineral_comp_names_.size();

  chem_engine_->GetAqueousKineticNames(aqueous_kinetics_names_);
  number_aqueous_kinetics_ = aqueous_kinetics_names_.size();

  chem_engine_->GetSurfaceSiteNames(sorption_site_names_);
  number_sorption_sites_ = sorption_site_names_.size();
  int num_sorption_components = chem_engine_->NumSorbedSpecies();
  // Amanzi only supports all or nothing sorption
  AMANZI_ASSERT(num_sorption_components == 0 ||
                num_sorption_components == number_aqueous_components_);

  chem_engine_->GetIonExchangeNames(ion_exchange_site_names_);
  number_ion_exchange_sites_ = ion_exchange_site_names_.size();

  int num_isotherm_components = chem_engine_->NumIsothermSpecies();
  // Amanzi only supports all or nothing isotherms
  AMANZI_ASSERT(num_isotherm_components == 0 ||
                num_isotherm_components == number_aqueous_components_);

  number_aux_data_ = chem_engine_->Sizes().num_aux_integers + chem_engine_->Sizes().num_aux_doubles;

  // things we always need
  dens_key_ = Keys::readKey(*plist_, domain_, "mass density liquid", "mass_density_liquid");
  poro_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
  sat_key_ = Keys::readKey(*plist_, domain_, "saturation liquid", "saturation_liquid");

  // NOTE: defaults to no temperature!  User can say to use temperature by setting this key.
  temp_key_ = Keys::readKey(*plist_, domain_, "temperature", "");

  if (number_mineral_components_ > 0) {
    mineral_volume_fraction_key_ = Keys::readKey(*plist_, domain_, "mineral volume fractions", "mineral_volume_fractions");
    mineral_specific_surface_area_key_ = Keys::readKey(*plist_, domain_, "mineral specific surface area", "mineral_specific_surface_area");
    mineral_rate_constant_key_ = Keys::readKey(*plist_, domain_, "mineral rate constant", "mineral_rate_constant");
  }

  if (number_sorption_sites_ > 0) {
    sorp_site_density_key_ = Keys::readKey(*plist_, domain_, "sorption site density", "sorption_site_density");
  }

  if (num_sorption_components > 0) {
    total_sorbed_key_ = Keys::readKey(*plist_, domain_, "total sorbed", "total_sorbed");
  }

  if (num_isotherm_components > 0) {
    // all of these are parameters
    isotherm_kd_key_ = Keys::readKey(*plist_, domain_, "isotherm kd", "isotherm_kd");
    isotherm_freundlich_n_key_ = Keys::readKey(*plist_, domain_, "isotherm freundlich n", "isotherm_freundlich_n");
    isotherm_langmuir_b_key_ = Keys::readKey(*plist_, domain_, "isotherm langmuir b", "isotherm_langmuir_b");
  }

  if (number_aqueous_kinetics_ > 0) {
    aqueous_kinetic_rate_constant_key_ = Keys::readKey(*plist_, domain_, "aqueous kinetic rate constant", "aqueous_kinetic_rate_constant");
  }

  if (number_ion_exchange_sites_ > 0) {
    cation_exchange_capacity_key_ = Keys::readKey(*plist_, domain_, "cation exchange capacity", "cation_exchange_capacity");
  }

  if (number_aux_data_ > 0) {
    aux_data_key_ = Keys::readKey(*plist_, domain_, "aux data", "aux_data");
  }

  // Set up auxiliary chemistry data using the ChemistryEngine.
  std::vector<std::string> aux_out_names;
  chem_engine_->GetAuxiliaryOutputNames(aux_out_names, aux_out_subfield_names_);

  for (const auto& name : aux_out_names) {
    if (name == "pH") {
      pH_key_ = Keys::readKey(*plist_, domain_, "pH", "pH");
      aux_out_names_[name] = pH_key_;
    } else if (name == "mineral_saturation_index") {
      mineral_sat_index_key_ = Keys::readKey(*plist_, domain_, "mineral saturation index", "mineral_saturation_index");
      aux_out_names_[name] = mineral_sat_index_key_;
    } else if (name == "mineral_reaction_rate") {
      mineral_reaction_rate_key_ = Keys::readKey(*plist_, domain_, "mineral reaction rate", "mineral_reaction_rate");
      aux_out_names_[name] = mineral_reaction_rate_key_;
    } else if (name == "aqueous_kinetic_rate") {
      aqueous_kinetic_rate_key_ = Keys::readKey(*plist_, domain_, "aqueous kinetic rate", "aqueous_kinetic_rate");
      aux_out_names_[name] = aqueous_kinetic_rate_key_;
    } else if (name == "primary_free_ion_concentration") {
      primary_ion_conc_key_ = Keys::readKey(*plist_, domain_, "primary free ion concentration", "primary_free_ion_concentration");
      aux_out_names_[name] = primary_ion_conc_key_;
    } else if (name == "primary_activity_coeff") {
      primary_activity_coef_key_ = Keys::readKey(*plist_, domain_, "primary activity coefficient", "primary_activity_coefficient");
      aux_out_names_[name] = primary_activity_coef_key_;
    } else if (name == "secondary_free_ion_concentration") {
      secondary_ion_conc_key_ = Keys::readKey(*plist_, domain_, "secondary free ion concentration", "secondary_free_ion_concentration");
      aux_out_names_[name] = secondary_ion_conc_key_;
    } else if (name == "secondary_activity_coeff") {
      secondary_activity_coef_key_ = Keys::readKey(*plist_, domain_, "secondary activity coefficient", "secondary_activity_coefficient");
      aux_out_names_[name] = secondary_activity_coef_key_;
    // } else if (name == "gas_partial_pressure") {
    //   gas_partial_pressure_key_ = Keys::readKey(*plist_, domain_, "gas partial pressure", "gas_partial_pressure");
    //   aux_out_names_[name] = gas_partial_pressure_key_;
    } else {
      // update this logic, keys, and CopyFromAlquimia(), adding new variables for storing diagnostics
      AMANZI_ASSERT(false);
    }
  }
}


/* ******************************************************************
* Register fields and evaluators with the State
******************************************************************* */
void
Alquimia_PK::Setup()
{
  Chemistry_PK::Setup();

  {  // external fields, only at the current time
    auto keys = std::vector<Key>{ dens_key_, poro_key_, sat_key_, temp_key_ };
    for (const auto& key : keys) {
      if (!key.empty()) {
        requireAtCurrent(key, tag_current_, *S_)
          .SetMesh(mesh_)
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      }
    }
  }

  { // parameters, only at the current
    auto keys = std::vector<Key>{ isotherm_kd_key_, isotherm_freundlich_n_key_, isotherm_langmuir_b_key_ };
    for (const auto& key : keys) {
      if (!key.empty()) {
        requireAtCurrent(key, tag_current_, *S_)
          .SetMesh(mesh_)
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_components_);
        S_->GetRecordSetW(key).set_subfieldnames(aqueous_comp_names_);
      }
    }

    if (!mineral_rate_constant_key_.empty()) {
      requireAtCurrent(mineral_rate_constant_key_, tag_current_, *S_)
        .SetMesh(mesh_)
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, number_mineral_components_);
      S_->GetRecordSetW(mineral_rate_constant_key_).set_subfieldnames(mineral_comp_names_);
    }

    if (!aqueous_kinetic_rate_constant_key_.empty()) {
      requireAtCurrent(aqueous_kinetic_rate_constant_key_, tag_current_, *S_)
        .SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_kinetics_);
      S_->GetRecordSetW(aqueous_kinetic_rate_constant_key_).set_subfieldnames(aqueous_kinetics_names_);
    }
  }

  { // aqueous state data
    auto keys = std::vector<Key>{key_, total_sorbed_key_};
    for (const auto& key : keys) {
      if (!key.empty()) {
        requireAtNext(key, tag_next_, *S_, true, passwd_)
          .SetMesh(mesh_)
          ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aqueous_components_);
        S_->GetRecordSetW(key).set_subfieldnames(aqueous_comp_names_);
        requireAtCurrent(key, tag_current_, *S_, passwd_);
      }
    }
  }

  { // mineral state data
    auto keys = std::vector<Key>{mineral_volume_fraction_key_, mineral_specific_surface_area_key_};
    for (const auto& key : keys) {
      if (!key.empty()) {
        requireAtNext(key, tag_next_, *S_, true, passwd_)
          .SetMesh(mesh_)
          ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_mineral_components_);
        S_->GetRecordSetW(key).set_subfieldnames(mineral_comp_names_);
        requireAtCurrent(key, tag_current_, *S_, passwd_);
      }
    }
  }

  // sorption site state data
  if (!sorp_site_density_key_.empty()) {
    requireAtNext(sorp_site_density_key_, tag_next_, *S_, true, passwd_)
      .SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_sorption_sites_);
    S_->GetRecordSetW(sorp_site_density_key_).set_subfieldnames(sorption_site_names_);
    requireAtCurrent(sorp_site_density_key_, tag_current_, *S_, passwd_);
  }

  // ion exchange site density
  if (!cation_exchange_capacity_key_.empty()) {
    requireAtNext(cation_exchange_capacity_key_, tag_next_, *S_, true, passwd_)
      .SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_ion_exchange_sites_);
    S_->GetRecordSetW(cation_exchange_capacity_key_).set_subfieldnames(ion_exchange_site_names_);
    requireAtCurrent(cation_exchange_capacity_key_, tag_current_, *S_, passwd_);
  }

  // aux_data
  if (!aux_data_key_.empty()){
    requireAtNext(aux_data_key_, tag_next_, *S_, true, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, number_aux_data_);
    requireAtCurrent(aux_data_key_, tag_current_, *S_, passwd_);
    S_->GetRecordW(aux_data_key_, tag_next_, passwd_).set_io_vis(false);
  }

  // aux output data are done one-at-time because they are different sizes.
  //
  // We only need these at the NEXT time as they are diagnostics, not state.
  //
  // No need to checkpoint these, despite being primary variables
  unsigned i = 0;
  for (const auto& [name, key] : aux_out_names_) {
    AMANZI_ASSERT(!key.empty());
    requireAtNext(key, tag_next_, *S_, true, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, aux_out_subfield_names_[i].size());
    S_->GetRecordW(key, tag_next_, passwd_).set_io_checkpoint(false);
    S_->GetRecordSetW(key).set_subfieldnames(aux_out_subfield_names_[i]);
    ++i;
  }
}


/* *******************************************************************
* Initialization
******************************************************************* */
void
Alquimia_PK::Initialize()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // initialization using the base class, will initialize tcc
  Chemistry_PK::Initialize();

  // Read XML parameters from our input file.
  XMLParameters();

  for (const auto& [name, key] : aux_out_names_) {
    S_->GetRecordW(key, tag_next_, passwd_).set_initialized();
  }

  if (std::abs(initial_conditions_time_ - S_->get_time()) < 1e-8 * (1.0 + fabs(S_->get_time()))) {
    updateSubstate();
    int ierr = 0;

    for (const auto& [region, condition] : chem_initial_conditions_) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        *vo_->os() << "enforcing geochemical condition \"" << condition << "\" in region \""
                   << region << "\"\n";
      }

      // Get the cells that belong to this region.
      int num_cells =
        mesh_->getSetSize(region, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
      auto cell_indices = mesh_->getSetEntities(
        region, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

      // Loop over the cells.
      for (int i = 0; i < num_cells; ++i) {
        int cell = cell_indices[i];
        ierr = initializeSingleCell_(cell, condition);
        if (ierr) break;
      }
    }

    // figure out if any of the processes threw an error, if so all processes will re-throw
    int recv = 0;
    mesh_->getComm()->MaxAll(&ierr, &recv, 1);
    if (recv != 0) {
      Errors::Message msg("Error in Alquimia_PK::Initialize()");
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* *******************************************************************
* This helper performs initialization on a single cell within Amanzi's state.
* It returns an error code that indicates success (0) or failure (1).
******************************************************************* */
int
Alquimia_PK::initializeSingleCell_(int cell, const std::string& condition)
{
  copyToAlquimia(cell, beaker_);
  chem_engine_->EnforceCondition(condition, S_->get_time(tag_current_),
          beaker_.properties, beaker_.state, beaker_.aux_data, beaker_.aux_output);
  copyFromAlquimia_(cell);

  // ETC: hacking to get consistent solution -- if there is no water
  // (e.g. surface system), we still need to call EnforceCondition() as it also
  // gets aux data set up correctly.  But the concentrations need to be
  // overwritten as 0 to get expected output.  Therefore we manually overwrite
  // this now.
  if (beaker_.properties.saturation <= saturation_tolerance_ && substate_.tcc_new) {
    for (int i = 0; i != substate_.tcc_new->NumVectors(); ++i) {
      (*substate_.tcc_new)[i][cell] = 0.;
    }
  }
  return 0;
}


/* *******************************************************************
 * This helper extracts a pointer to Epetra_MultiVector into the substate
 * struct for use in AdvanceSingleCell.
 ******************************************************************* */
void
Alquimia_PK::updateSubstate()
{
  // external things
  S_->GetEvaluator(poro_key_, tag_current_).Update(*S_, name_);
  substate_.porosity = &*S_->Get<CompositeVector>(poro_key_, tag_current_).ViewComponent("cell", false);

  S_->GetEvaluator(sat_key_, tag_current_).Update(*S_, name_);
  substate_.saturation_liquid = &*S_->Get<CompositeVector>(sat_key_, tag_current_).ViewComponent("cell", false);

  S_->GetEvaluator(dens_key_, tag_current_).Update(*S_, name_);
  substate_.mass_density_liquid = &*S_->Get<CompositeVector>(dens_key_, tag_current_).ViewComponent("cell", false);

  if (!temp_key_.empty()) {
    S_->GetEvaluator(temp_key_, tag_current_).Update(*S_, name_);
    substate_.temperature = &*S_->Get<CompositeVector>(temp_key_, tag_current_).ViewComponent("cell", false);
  }

  // parameters
  if (!isotherm_kd_key_.empty()) {
    S_->GetEvaluator(isotherm_kd_key_, tag_current_).Update(*S_, name_);
    substate_.isotherm_kd = &*S_->Get<CompositeVector>(isotherm_kd_key_, tag_current_).ViewComponent("cell", false);

    S_->GetEvaluator(isotherm_freundlich_n_key_, tag_current_).Update(*S_, name_);
    substate_.isotherm_freundlich_n = &*S_->Get<CompositeVector>(isotherm_freundlich_n_key_, tag_current_).ViewComponent("cell", false);

    S_->GetEvaluator(isotherm_langmuir_b_key_, tag_current_).Update(*S_, name_);
    substate_.isotherm_langmuir_b = &*S_->Get<CompositeVector>(isotherm_langmuir_b_key_, tag_current_).ViewComponent("cell", false);
  }

  if (!mineral_rate_constant_key_.empty()) {
    S_->GetEvaluator(mineral_rate_constant_key_, tag_current_).Update(*S_, name_);
    substate_.mineral_rate_constant = &*S_->Get<CompositeVector>(mineral_rate_constant_key_, tag_current_).ViewComponent("cell", false);
  }

  if (!aqueous_kinetic_rate_constant_key_.empty()) {
    S_->GetEvaluator(aqueous_kinetic_rate_constant_key_, tag_current_).Update(*S_, name_);
    substate_.aqueous_kinetic_rate_constant = &*S_->Get<CompositeVector>(aqueous_kinetic_rate_constant_key_, tag_current_).ViewComponent("cell", false);
  }

  // internal things
  if (operator_split_) {
    // OPERATOR SPLITTING -- TCC old is NEXT here!
    substate_.tcc_old = &*S_->Get<CompositeVector>(key_, tag_next_).ViewComponent("cell", false);
  } else {
    // stand-alone chemistry -- TCC is CURRENT here!
    substate_.tcc_old = &*S_->Get<CompositeVector>(key_, tag_current_).ViewComponent("cell", false);
  }
  substate_.tcc_new = &*S_->GetW<CompositeVector>(key_, tag_next_, passwd_).ViewComponent("cell", false);

  if (number_mineral_components_ > 0) {
    substate_.mineral_volume_fraction_old = &*S_->Get<CompositeVector>(mineral_volume_fraction_key_, tag_current_).ViewComponent("cell", false);
    substate_.mineral_volume_fraction_new = &*S_->GetW<CompositeVector>(mineral_volume_fraction_key_, tag_next_, passwd_).ViewComponent("cell", false);

    substate_.mineral_specific_surface_area_old = &*S_->Get<CompositeVector>(mineral_specific_surface_area_key_, tag_current_).ViewComponent("cell", false);
    substate_.mineral_specific_surface_area_new = &*S_->GetW<CompositeVector>(mineral_specific_surface_area_key_, tag_next_, passwd_).ViewComponent("cell", false);
  }

  if (!total_sorbed_key_.empty()) {
    substate_.total_sorbed_old = &*S_->Get<CompositeVector>(total_sorbed_key_, tag_current_).ViewComponent("cell", false);
    substate_.total_sorbed_new = &*S_->GetW<CompositeVector>(total_sorbed_key_, tag_next_, passwd_).ViewComponent("cell", false);
  }

  if (number_ion_exchange_sites_ > 0) {
    substate_.cation_exchange_capacity_old = &*S_->Get<CompositeVector>(cation_exchange_capacity_key_, tag_current_).ViewComponent("cell", false);
    substate_.cation_exchange_capacity_old = &*S_->GetW<CompositeVector>(cation_exchange_capacity_key_, tag_next_, passwd_).ViewComponent("cell", false);
  }

  // aux data
  if (!aux_data_key_.empty()) {
    substate_.aux_data_old = &*S_->Get<CompositeVector>(aux_data_key_, tag_current_).ViewComponent("cell", false);
    substate_.aux_data_new = &*S_->GetW<CompositeVector>(aux_data_key_, tag_next_, passwd_).ViewComponent("cell", false);
  }

  // aux output
  if (!pH_key_.empty())
    substate_.pH = &*S_->GetW<CompositeVector>(pH_key_, tag_next_, passwd_).ViewComponent("cell", false);
  if (!mineral_sat_index_key_.empty())
    substate_.mineral_saturation_index = &*S_->GetW<CompositeVector>(mineral_sat_index_key_, tag_next_, passwd_).ViewComponent("cell", false);
  if (!mineral_reaction_rate_key_.empty())
    substate_.mineral_reaction_rate = &*S_->GetW<CompositeVector>(mineral_reaction_rate_key_, tag_next_, passwd_).ViewComponent("cell", false);
  if (!aqueous_kinetic_rate_key_.empty())
    substate_.aqueous_kinetic_rate = &*S_->GetW<CompositeVector>(aqueous_kinetic_rate_key_, tag_next_, passwd_).ViewComponent("cell", false);
  if (!primary_ion_conc_key_.empty())
    substate_.primary_free_ion_concentration = &*S_->GetW<CompositeVector>(primary_ion_conc_key_, tag_next_, passwd_).ViewComponent("cell", false);
  if (!primary_activity_coef_key_.empty())
    substate_.primary_activity_coefficient = &*S_->GetW<CompositeVector>(primary_activity_coef_key_, tag_next_, passwd_).ViewComponent("cell", false);
  if (!secondary_ion_conc_key_.empty())
    substate_.secondary_free_ion_concentration = &*S_->GetW<CompositeVector>(secondary_ion_conc_key_, tag_next_, passwd_).ViewComponent("cell", false);
  if (!secondary_activity_coef_key_.empty())
    substate_.secondary_activity_coefficient = &*S_->GetW<CompositeVector>(secondary_activity_coef_key_, tag_next_, passwd_).ViewComponent("cell", false);
  // if (!gas_partial_pressure_key_.empty())
  //   substate_.gas_partial_pressure = &*S_->GetW<CompositeVector>(gas_partial_pressure_key_, tag_next_, passwd_).ViewComponent("cell", false);

}


/* *******************************************************************
* Initialization helper functions
******************************************************************* */
void
Alquimia_PK::ParseChemicalConditionRegions_(const Teuchos::ParameterList& param_list,
                                           std::map<std::string, std::string>& conditions)
{
  Errors::Message msg;

  // Go through the sublist containing the chemical conditions.
  for (auto it = param_list.begin(); it != param_list.end(); ++it) {
    // This parameter list contains sublists, each corresponding to a condition.
    std::string cond_name = param_list.name(it);
    assert(param_list.isSublist(cond_name));
    const Teuchos::ParameterList& cond_sublist = param_list.sublist(cond_name);

    // Apply this condition to all desired regions.
    if (!cond_sublist.isType<Teuchos::Array<std::string>>("regions")) {
      msg << "Alquimia_PK::ParseChemicalConditionRegions(): \n";
      msg << "  Geochemical condition '" << cond_name << "' has no valid 'regions' entry.\n";
      Exceptions::amanzi_throw(msg);
    }

    Teuchos::Array<std::string> regions = cond_sublist.get<Teuchos::Array<std::string>>("regions");
    for (size_t r = 0; r < regions.size(); ++r) {
      // We allow for cell-based and face-based regions to accommodate both
      // initial and boundary conditions.
      if (!mesh_->isValidSetName(regions[r], AmanziMesh::Entity_kind::CELL) &&
          !mesh_->isValidSetName(regions[r], AmanziMesh::Entity_kind::FACE)) {
        msg << "Alquimia_PK::ParseChemicalConditionRegions(): \n";
        msg << "  Invalid region '" << regions[r] << "' given for geochemical condition '"
            << cond_name << "'.\n";
        Exceptions::amanzi_throw(msg);
      }
      conditions[regions[r]] = cond_name;
    }
  }
}


/* *******************************************************************
*
******************************************************************* */
void
Alquimia_PK::XMLParameters()
{
  Errors::Message msg;
  Teuchos::OSTab tab = vo_->getOSTab();

  // Add any geochemical conditions we find in the Chemistry section of the file.
  if (plist_->isParameter("geochemical conditions")) {
    Teuchos::ParameterList conditions = plist_->sublist("geochemical conditions");
    for (auto it = conditions.begin(); it != conditions.end(); ++it) {
      // This parameter list contains sublists, each corresponding to a
      // geochemical condition.
      std::string cond_name = conditions.name(it);
      assert(conditions.isSublist(cond_name));
      const Teuchos::ParameterList& cond_sublist = conditions.sublist(cond_name);

      // Create the entry for this geochemical condition within the chemistry engine,
      // overwriting any previous definition.
      chem_engine_->CreateCondition(cond_name);

      // Now mine the entry for details.
      for (auto it2 = cond_sublist.begin(); it2 != cond_sublist.end(); ++it2) {
        std::string species_name = cond_sublist.name(it2);
        assert(cond_sublist.isSublist(species_name));
        const Teuchos::ParameterList& aqueous_constraint = cond_sublist.sublist(species_name);

        // If the primary species has an associated equilibration constraint, we need to retrieve its information.
        std::string equilibrate_name;
        if (aqueous_constraint.isParameter("equilibrate")) {
          equilibrate_name = aqueous_constraint.get<std::string>("equilibrate");

          // The information for this mineral species must appear alongside the
          // aqueous constraint in the geochemical condition.
        }

        // What kind of aqueous constraint do we have on this species?
        static const char* valid_types[] = { "total_aqueous", "total_sorb", "free",  "mineral",
                                             "gas",           "pH",         "charge" };
        static int num_valid_types = 7;
        std::string type;
        for (int i = 0; i < num_valid_types; ++i) {
          if (aqueous_constraint.isParameter(valid_types[i])) type = std::string(valid_types[i]);
          // check whether a mineral or gas name have been provided to equilibrate with
          if (type == "mineral" || type == "gas") {
            if (equilibrate_name == "") {
              msg << "No mineral or gas species has been given to equilibrate with species '"
                  << species_name << "'.\n";
              Exceptions::amanzi_throw(msg);
            }
          }
        }
        if (!type.empty()) {
          // It's a valid aqueous constraint, so we add it to the chemistry engine
          // under the current geochemical condition.
          chem_engine_->AddAqueousConstraint(cond_name, species_name, type, equilibrate_name);
        } else {
          // We have an invalid aqueous contraint.
          msg << "Invalid aqueous constraint type for " << species_name << ".\n";
          msg << "Valid types are total_aqueous, total_sorb, free, mineral, gas, pH, and charge.\n";
          Exceptions::amanzi_throw(msg);
        }
      }
    }
  }

  // Now associate regions with chemical conditions based on initial
  // condition specifications in the file.
  auto initial_conditions = Teuchos::sublist(plist_, "initial conditions");
  if (initial_conditions->isSublist("geochemical conditions")) {
    Teuchos::ParameterList& geochem_conditions =
      initial_conditions->sublist("geochemical conditions");
    ParseChemicalConditionRegions_(geochem_conditions, chem_initial_conditions_);
    if (chem_initial_conditions_.empty()) {
      if (plist_->isSublist("initial conditions")) {
        msg << "Alquimia_PK::XMLParameters(): No geochemical conditions were found in "
               "\"PK->initial condition->geochemical conditions\"";
      } else {
        msg << "Alquimia_PK::XMLParameters(): No geochemical conditions were found in "
               "\"State->initial conditions->geochemical conditions\"";
      }
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* *******************************************************************
*
******************************************************************* */
void
Alquimia_PK::copyToAlquimia(int cell, AlquimiaBeaker& beaker)
{
  beaker.properties.volume = mesh_->getCellVolume(cell);
  beaker.properties.saturation = (*substate_.saturation_liquid)[0][cell];
  beaker.state.water_density = (*substate_.mass_density_liquid)[0][cell];
  beaker.state.porosity = (*substate_.porosity)[0][cell];

  if (substate_.temperature) {
    beaker.state.temperature = (*substate_.temperature)[0][cell];
  }

  for (unsigned int i = 0; i != number_aqueous_components_; ++i) {
    beaker.state.total_mobile.data[i] = (*substate_.tcc_old)[i][cell];

    if (substate_.total_sorbed_old) {
      beaker.state.total_immobile.data[i] = (*substate_.total_sorbed_old)[i][cell];
    }
  }

  // minerals
  for (unsigned int i = 0; i != number_mineral_components_; ++i) {
    beaker.state.mineral_volume_fraction.data[i] = (*substate_.mineral_volume_fraction_old)[i][cell];
    beaker.properties.mineral_rate_cnst.data[i] = (*substate_.mineral_rate_constant)[i][cell];
    beaker.state.mineral_specific_surface_area.data[i] = (*substate_.mineral_specific_surface_area_old)[i][cell];
  }

  // surface complexation
  for (unsigned int i = 0; i != number_sorption_sites_; ++i)
    beaker.state.surface_site_density.data[i] = (*substate_.sorption_site_density_old)[i][cell];

  // ion exchange
  for (unsigned int i = 0; i != number_ion_exchange_sites_; ++i)
    beaker.state.cation_exchange_capacity.data[i] = (*substate_.cation_exchange_capacity_old)[i][cell];

  // sorption isotherms
  if (substate_.isotherm_kd) {
    for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
      beaker.properties.isotherm_kd.data[i] = (*substate_.isotherm_kd)[i][cell];
      beaker.properties.freundlich_n.data[i] = (*substate_.isotherm_freundlich_n)[i][cell];
      beaker.properties.langmuir_b.data[i] = (*substate_.isotherm_langmuir_b)[i][cell];
    }
  }

  // aqueous kinetics
  for (unsigned int i = 0; i != number_aqueous_kinetics_; ++i) {
    beaker.properties.aqueous_kinetic_rate_cnst.data[i] = (*substate_.aqueous_kinetic_rate_constant)[i][cell];
  }

  // Auxiliary data
  int num_aux_ints = chem_engine_->Sizes().num_aux_integers;
  int num_aux_doubles = chem_engine_->Sizes().num_aux_doubles;
  for (unsigned int i = 0; i != num_aux_ints; ++i) {
    beaker.aux_data.aux_ints.data[i] = (int)(*substate_.aux_data_old)[i][cell];
  }
  for (unsigned int i = 0; i != num_aux_doubles; ++i) {
    beaker.aux_data.aux_doubles.data[i] = (*substate_.aux_data_old)[i+num_aux_ints][cell];
  }
}


/* *******************************************************************
*
******************************************************************* */
void
Alquimia_PK::copyFromAlquimia_(int cell)
{
  for (unsigned i = 0; i != number_aqueous_components_; ++i) {
    (*substate_.tcc_new)[i][cell] = beaker_.state.total_mobile.data[i];

    if (substate_.total_sorbed_new) {
      (*substate_.total_sorbed_new)[i][cell] = beaker_.state.total_immobile.data[i];
    }
  }

  // Mineral properties.
  for (unsigned i = 0; i != number_mineral_components_; ++i) {
    (*substate_.mineral_volume_fraction_new)[i][cell] = beaker_.state.mineral_volume_fraction.data[i];
    (*substate_.mineral_specific_surface_area_new)[i][cell] = beaker_.state.mineral_specific_surface_area.data[i];
  }

  // surface complexation
  for (unsigned i = 0; i != number_sorption_sites_; ++i) {
    (*substate_.sorption_site_density_new)[i][cell] = beaker_.state.surface_site_density.data[i];
  }

  // ion exchange
  for (unsigned i = 0; i != number_ion_exchange_sites_; ++i) {
    (*substate_.cation_exchange_capacity_new)[i][cell] = beaker_.state.cation_exchange_capacity.data[i];
  }

  int num_aux_ints = chem_engine_->Sizes().num_aux_integers;
  int num_aux_doubles = chem_engine_->Sizes().num_aux_doubles;
  for (unsigned i = 0; i != num_aux_ints; ++i) {
    (*substate_.aux_data_new)[i][cell] = (double)beaker_.aux_data.aux_ints.data[i];
  }
  for (unsigned i = 0; i != num_aux_doubles; ++i) {
    (*substate_.aux_data_new)[i+num_aux_ints][cell] = beaker_.aux_data.aux_doubles.data[i];
  }

  if (substate_.pH)
    (*substate_.pH)[0][cell] = beaker_.aux_output.pH;

  if (substate_.mineral_saturation_index)
    for (unsigned i = 0; i != substate_.mineral_saturation_index->NumVectors(); ++i)
      (*substate_.mineral_saturation_index)[i][cell] = beaker_.aux_output.mineral_saturation_index.data[i];

  if (substate_.mineral_reaction_rate)
    for (unsigned i = 0; i != substate_.mineral_reaction_rate->NumVectors(); ++i)
      (*substate_.mineral_reaction_rate)[i][cell] = beaker_.aux_output.mineral_reaction_rate.data[i];

  if (substate_.aqueous_kinetic_rate)
    for (unsigned i = 0; i != substate_.aqueous_kinetic_rate->NumVectors(); ++i)
      (*substate_.aqueous_kinetic_rate)[i][cell] = beaker_.aux_output.aqueous_kinetic_rate.data[i];

  if (substate_.primary_free_ion_concentration)
    for (unsigned i = 0; i != substate_.primary_free_ion_concentration->NumVectors(); ++i)
      (*substate_.primary_free_ion_concentration)[i][cell] = beaker_.aux_output.primary_free_ion_concentration.data[i];

  if (substate_.primary_activity_coefficient)
    for (unsigned i = 0; i != substate_.primary_activity_coefficient->NumVectors(); ++i)
      (*substate_.primary_activity_coefficient)[i][cell] = beaker_.aux_output.primary_activity_coeff.data[i];

  if (substate_.secondary_free_ion_concentration)
    for (unsigned i = 0; i != substate_.secondary_free_ion_concentration->NumVectors(); ++i)
      (*substate_.secondary_free_ion_concentration)[i][cell] = beaker_.aux_output.secondary_free_ion_concentration.data[i];

  if (substate_.secondary_activity_coefficient)
    for (unsigned i = 0; i != substate_.secondary_activity_coefficient->NumVectors(); ++i)
      (*substate_.secondary_activity_coefficient)[i][cell] = beaker_.aux_output.secondary_activity_coeff.data[i];

  // need num gas, which isn't yet used... --ETC
  // if (substate_.gas_partial_pressure)
  //   (*substate_.gas_partial_pressure)[0][cell] = beaker_.aux_output.gas_partial_pressure;
}


/* *******************************************************************
* This helper advances the solution on a single cell within Amanzi's state.
* It returns the number of iterations taken to obtain the advanced solution,
* or -1 if an error occurred.
******************************************************************* */
int
Alquimia_PK::advanceSingleCell_(int cell, double dt)
{
  // Copy the state and property information from Amanzi's state within
  // this cell to Alquimia.
  copyToAlquimia(cell, beaker_);

  int num_iterations = 0;
  if (beaker_.properties.saturation > saturation_tolerance_) {
    bool success =
      chem_engine_->Advance(dt,
                            beaker_.properties,
                            beaker_.state,
                            beaker_.aux_data,
                            beaker_.aux_output,
                            num_iterations,
                            mesh_->getMap(AmanziMesh::Entity_kind::CELL, false).GID(cell));
    if (not success) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "no convergence in cell: "
                   << mesh_->getMap(AmanziMesh::Entity_kind::CELL, false).GID(cell) << std::endl;
      }
      return -1;
    }

    // Copy the information back into Amanzi's state, updating the given total
    // concentration vector.
    copyFromAlquimia_(cell);
  }
  return num_iterations;
}


/* *******************************************************************
* Copies all state data from tag_source to tag_dest.
******************************************************************* */
void
Alquimia_PK::copyFields_(const Tag& tag_dest, const Tag& tag_source)
{

  auto keys = std::vector<Key>{ mineral_volume_fraction_key_,
    mineral_specific_surface_area_key_,
    sorp_site_density_key_,
    cation_exchange_capacity_key_,
    aux_data_key_ };

  for (const auto& key : keys) {
    if (!key.empty()) assign(key, tag_dest, tag_source, *S_);
  }
}



} // namespace AmanziChemistry
} // namespace Amanzi
