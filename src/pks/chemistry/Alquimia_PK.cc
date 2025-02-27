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

// Chemistry
#include "Alquimia_PK.hh"

namespace Amanzi {
namespace AmanziChemistry {


/* *******************************************************************
* Constructor
******************************************************************* */
Alquimia_PK::Alquimia_PK(Teuchos::ParameterList& pk_tree,
                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, glist, S, soln),
    Chemistry_PK(pk_tree, glist, S, soln),
    chem_initialized_(false)
{}


/* *******************************************************************
* Destroy ansilary data structures.
******************************************************************* */
Alquimia_PK::~Alquimia_PK()
{
  if (chem_initialized_)
    chem_engine_->FreeState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);
}


/* *******************************************************************
* Parser
******************************************************************* */
void
Alquimia_PK::parseParameterList()
{
  Chemistry_PK::parseParameterList();

  // create chemistry engine. (should we do it later in Setup()?)
  if (!plist_->isParameter("engine")) {
    Errors::Message msg;
    msg << "No 'engine' parameter found in the parameter list for 'Chemistry'.\n";
    Exceptions::amanzi_throw(msg);
  }
  if (!plist_->isParameter("engine input file")) {
    Errors::Message msg;
    msg << "No 'engine input file' parameter found in the parameter list for 'Chemistry'.\n";
    Exceptions::amanzi_throw(msg);
  }
  std::string engine_name = plist_->get<std::string>("engine");
  std::string engine_inputfile = plist_->get<std::string>("engine input file");
  chem_engine_ = Teuchos::rcp(new AmanziChemistry::ChemistryEngine(engine_name, engine_inputfile));

  // grab the component names
  comp_names_.clear();
  chem_engine_->GetPrimarySpeciesNames(comp_names_);

  number_aqueous_components_ = comp_names_.size();
  number_free_ion_ = number_aqueous_components_;
  number_total_sorbed_ = number_aqueous_components_;

  chem_engine_->GetAqueousKineticNames(aqueous_kinetics_names_);
  number_aqueous_kinetics_ = aqueous_kinetics_names_.size();
}


/* ******************************************************************
* Register fields and evaluators with the State
******************************************************************* */
void
Alquimia_PK::Setup()
{
  Chemistry_PK::Setup();

  // Set up auxiliary chemistry data using the ChemistryEngine.
  chem_engine_->GetAuxiliaryOutputNames(aux_names_, aux_subfield_names_);

  for (size_t i = 0; i < aux_names_.size(); ++i) {
    aux_names_[i] = Keys::getKey(domain_, aux_names_[i]);

    if (!S_->HasRecord(aux_names_[i])) {
      S_->Require<CompositeVector, CompositeVectorSpace>(
          aux_names_[i], tag_next_, passwd_, aux_subfield_names_[i])
        .SetMesh(mesh_)
        ->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, aux_subfield_names_[i].size());
    }
  }

  if (plist_->isParameter("auxiliary data")) {
    auto names = plist_->get<Teuchos::Array<std::string>>("auxiliary data");

    for (auto it = names.begin(); it != names.end(); ++it) {
      Key aux_field_name = Keys::getKey(domain_, *it);
      if (!S_->HasRecord(aux_field_name)) {
        S_->Require<CompositeVector, CompositeVectorSpace>(aux_field_name, tag_next_, passwd_)
          .SetMesh(mesh_)
          ->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      }
    }
  }

  // Setup more auxiliary data
  if (!S_->HasRecord(alquimia_aux_data_key_, tag_next_)) {
    int num_aux_data =
      chem_engine_->Sizes().num_aux_integers + chem_engine_->Sizes().num_aux_doubles;
    S_->Require<CompositeVector, CompositeVectorSpace>(alquimia_aux_data_key_, tag_next_, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, num_aux_data);

    S_->GetRecordW(alquimia_aux_data_key_, tag_next_, passwd_).set_io_vis(false);
  }
}


/* *******************************************************************
* Initialization
******************************************************************* */
void
Alquimia_PK::Initialize()
{
  // initialization using the base class
  Chemistry_PK::Initialize();

  if (!aux_names_.empty()) {
    int n_total = 0;
    for (const auto& subfield_name : aux_subfield_names_) n_total += subfield_name.size();
    aux_output_ = Teuchos::rcp(
      new Epetra_MultiVector(mesh_->getMap(AmanziMesh::Entity_kind::CELL, false), n_total));
  } else {
    aux_output_ = Teuchos::null;
  }

  // Read XML parameters from our input file.
  XMLParameters();

  // initialize fields as soon as possible
  for (size_t i = 0; i < aux_names_.size(); ++i) {
    InitializeCVField(S_, *vo_, aux_names_[i], tag_next_, passwd_, 0.0);
  }

  // Initialize the data structures that we will use to traffic data between
  // Amanzi and Alquimia.
  chem_engine_->InitState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);

  if (using_sorption_ && alq_state_.total_immobile.data == NULL) {
    Errors::Message msg("Alquimia's state has no memory for total_immobile.");
    Exceptions::amanzi_throw(msg);
  }

  // all memory allocation consistency checks should be placed here
  AMANZI_ASSERT(alq_state_.surface_site_density.size == number_sorption_sites_);

  chem_engine_->GetMineralNames(mineral_names_);
  chem_engine_->GetPrimarySpeciesNames(primary_names_);
  InitializeAuxNamesMap_();

  // Do we need to initialize chemistry?
  int ierr = 0;

  // Ensure dependencies are filled
  S_->GetEvaluator(poro_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(fluid_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_key_, Tags::DEFAULT).Update(*S_, name_);

  if (fabs(initial_conditions_time_ - S_->get_time()) < 1e-8 * (1.0 + fabs(S_->get_time()))) {
    for (auto it = chem_initial_conditions_.begin(); it != chem_initial_conditions_.end(); ++it) {
      std::string region = it->first;
      std::string condition = it->second;

      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        Teuchos::OSTab tab = vo_->getOSTab();
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
        ierr = InitializeSingleCell(cell, condition);
      }
    }
  }

  // figure out if any of the processes threw an error, if so all processes will re-throw
  int recv = 0;
  mesh_->getComm()->MaxAll(&ierr, &recv, 1);
  if (recv != 0) {
    Errors::Message msg("Error in Alquimia_PK::Initialize()");
    Exceptions::amanzi_throw(msg);
  }

  // now publish auxiliary data to state
  if (aux_output_ != Teuchos::null) {
    int counter = 0;
    for (int i = 0; i < aux_names_.size(); ++i) {
      auto& aux_state =
        *S_->GetW<CompositeVector>(aux_names_[i], tag_next_, passwd_).ViewComponent("cell");
      for (int j = 0; j < aux_subfield_names_[i].size(); ++j) {
        *aux_state(j) = *(*aux_output_)(counter++);
      }
    }
  }

  chem_initialized_ = true;
}


/* *******************************************************************
* This helper performs initialization on a single cell within Amanzi's state.
* It returns an error code that indicates success (0) or failure (1).
******************************************************************* */
int
Alquimia_PK::InitializeSingleCell_(int cell, const std::string& condition)
{
  // NOTE: this should get set not to be hard-coded to Tags::DEFAULT, but
  // should use the same tag as transport.  See #673
  CopyToAlquimia(
    cell, aqueous_components_, alq_mat_props_, alq_state_, alq_aux_data_, Tags::DEFAULT);

  chem_engine_->EnforceCondition(
    condition, current_time_, alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);

  CopyAlquimiaStateToAmanzi(
    cell, alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_, aqueous_components_);


  // ETC: hacking to get consistent solution -- if there is no water
  // (e.g. surface system, we still need to call EnforceCondition() as it also
  // gets aux data set up correctly.  But the concentrations need to be
  // overwritten as 0 to get expected output.  Therefore we manually overwrite
  // this now.  Previously this happened due to a bug in ATS's reactive
  // transport coupler -- happy accidents.
  if (alq_mat_props_.saturation <= saturation_tolerance_)
    for (int i = 0; i != aqueous_components_->NumVectors(); ++i)
      (*aqueous_components_)[i][cell] = 0.;
  return 0;
}


Impl::AlquimiaSubstate
Alquimia_PK::createSubState()
{
  // external things
  porosity = &*S_->Get<CompositeVector>(poro_key_, tag_current_).ViewComponent("cell", false);
  water_saturation = &*S_->Get<CompositeVector>(saturation_key_, tag_current_).ViewComponent("cell", false);
  fluid_density = &*S_->Get<CompositeVector>(fluid_density_key_, tag_current_).ViewComponent("cell", false);

  if (S_->HasRecord(temperature_key_, tag_current_)) {
    temperature = &*S_->Get<CompositeVector>(temperature_key_, tag_current_).ViewComponent("cell", false);
  }

  // internal things
  if (number_aqueous_components_ > 0) {
    aqueous_components_old = &*S_->Get<CompositeVector>(key_, tag_current_).ViewComponent("cell", false);
    aqueous_components_new = &*S_->GetW<CompositeVector>(key_, tag_next_).ViewComponent("cell", false);

    free_ion_species_old = &*S_->Get<CompositeVector>(free_ion_species_key_, tag_current_).ViewComponent("cell", false);
    free_ion_species_new = &*S_->GetW<CompositeVector>(free_ion_species_key_, tag_next_).ViewComponent("cell", false);
  }

  if (number_mineral_components_ > 0) {
    mineral_vol_frac_old = &*S_->Get<CompositeVector>(mineral_vol_frac_key_, tag_current_).ViewComponent("cell", false);
    mineral_vol_frac_new = &*S_->GetW<CompositeVector>(mineral_vol_frac_key_, tag_next_).ViewComponent("cell", false);

    mineral_ssa_old = &*S_->Get<CompositeVector>(mineral_ssa_key_, tag_current_).ViewComponent("cell", false);
    mineral_ssa_new = &*S_->GetW<CompositeVector>(mineral_ssa_key_, tag_next_).ViewComponent("cell", false);

    mineral_rate_constant_old = &*S_->Get<CompositeVector>(mineral_rate_constant_key_, tag_current_).ViewComponent("cell", false);
    mineral_rate_constant_new = &*S_->GetW<CompositeVector>(mineral_rate_constant_key_, tag_next_).ViewComponent("cell", false);
  }

  if (using_sorption_) {
    total_sorbed_old = &*S_->Get<CompositeVector>(total_sorbed_key_, tag_current_).ViewComponent("cell", false);
    total_sorbed_new = &*S_->GetW<CompositeVector>(total_sorbed_key_, tag_next_).ViewComponent("cell", false);
  }

  if (number_ion_exchange_sites_ > 0) {
    ion_exchange_sites_old = &*S_->Get<CompositeVector>(ion_exchange_sites_key_, tag_current_).ViewComponent("cell", false);
    ion_exchange_sites_new = &*S_->GetW<CompositeVector>(ion_exchange_sites_key_, tag_next_).ViewComponent("cell", false);
  }

  if (using_sorption_isotherms_) {
    isotherm_kd_old = &*S_->Get<CompositeVector>(isotherm_kd_key_, tag_current_).ViewComponent("cell", false);
    isotherm_kd_new = &*S_->GetW<CompositeVector>(isotherm_kd_key_, tag_next_).ViewComponent("cell", false);

    isotherm_freundlich_n_old = &*S_->Get<CompositeVector>(isotherm_freundlich_n_key_, tag_current_).ViewComponent("cell", false);
    isotherm_freundlich_n_new = &*S_->GetW<CompositeVector>(isotherm_freundlich_n_key_, tag_next_).ViewComponent("cell", false);

    isotherm_langmuir_b_old = &*S_->Get<CompositeVector>(isotherm_langmuir_b_key_, tag_current_).ViewComponent("cell", false);
    isotherm_langmuir_b_new = &*S_->GetW<CompositeVector>(isotherm_langmuir_b_key_, tag_next_).ViewComponent("cell", false);
  }

  if (number_aqueous_kinetics_ > 0) {
    first_order_decay_constant_old = &*S_->Get<CompositeVector>(first_order_decay_constant_key_, tag_current_).ViewComponent("cell", false);
    first_order_decay_constant_new = &*S_->GetW<CompositeVector>(first_order_decay_constant_key_, tag_next_).ViewComponent("cell", false);
  }

  // aux data
  int num_aux_data = chem_engine_->Sizes().num_aux_integers + chem_engine_->Sizes().num_aux_doubles;
  if (num_aux_data > 0) {
    aux_data_old = &*S_->Get<CompositeVector>(aux_data_key_, tag_current_).ViewComponent("cell", false);
    aux_data_new = &*S_->GetW<CompositeVector>(aux_data_key_, tag_next_).ViewComponent("cell", false);
  }

  // aux output
  if (aux_subfield_names_.size() > 0) {
    aux_output = &*S_->GetW<CompositeVector>(aux_output_key_, tag_next_).ViewComponent("cell", false);
  }
}


/* *******************************************************************
* Initialization helper functions
******************************************************************* */
void
Alquimia_PK::ParseChemicalConditionRegions(const Teuchos::ParameterList& param_list,
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
  if (!glist_->isSublist("state")) {
    msg << "Alquimia_PK::XMLParameters(): \n";
    msg << "  No 'State' sublist was found!\n";
    Exceptions::amanzi_throw(msg);
  }

  Teuchos::RCP<Teuchos::ParameterList> initial_conditions;
  if (plist_->isSublist("initial condition")) {
    // ATS-style input spec -- initial conditions in the PK
    initial_conditions = Teuchos::sublist(plist_, "initial condition");
  } else {
    // Amanzi-style input spec -- initial conditions in State
    initial_conditions = Teuchos::sublist(Teuchos::sublist(glist_, "state"), "initial conditions");
  }
  if (initial_conditions->isSublist("geochemical conditions")) {
    Teuchos::ParameterList& geochem_conditions =
      initial_conditions->sublist("geochemical conditions");
    ParseChemicalConditionRegions(geochem_conditions, chem_initial_conditions_);
    if (chem_initial_conditions_.empty()) {
      if (plist_->isSublist("initial condition")) {
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
Alquimia_PK::CopyToAlquimia(const AlquimiaSubstate& state,
                            int cell)
{
  alq_mat_props_.volume = mesh_->getCellVolume(cell);
  alq_mat_props_.saturation = water_saturation[0][cell];
  alq_state_.water_density = fluid_density[0][cell];
  alq_state_.porosity = porosity[0][cell];

  if (state.temperature) {
    alq_state_.temperature = state.temperature[0][cell];
  }

  for (unsigned int i = 0; i != number_aqueous_components_; ++i) {
    alq_state_.total_mobile.data[i] = state.aqueous_components_old_old[i][cell];

    if (state.total_sorbed) {
      alq_state_.total_immobile.data[i] = state.total_sorbed_old[i][cell];
    }
  }

  // minerals
  for (unsigned int i = 0; i != number_minerals_; ++i) {
    alq_state_.mineral_volume_fraction.data[i] = state.mineral_volume_fraction_old[i][cell];
    alq_mat_props_.mineral_rate_cnst.data[i] = state.mineral_rate_constant_old[i][cell];
    alq_state_.mineral_specific_surface_area.data[i] = state.mineral_ssa_old[i][cell];
  }

  // ion exchange
  for (unsigned int i = 0; i != number_ion_exchange_sites_; ++i) {
    alq_state_.cation_exchange_capacity.data[i] = state.ion_exchange_old[i][cell];
  }

  // surface complexation
  for (unsigned int i = 0; i != number_sorption_sites_; ++i) {
    alq_state_.surface_site_density.data[i] = state.sorption_sites_old[i][cell];
  }

  // sorption isotherms
  if (using_sorption_isotherms_) {
    for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
      alq_mat_props_.isotherm_kd.data[i] = state.isotherm_kd_old[i][cell];
      alq_mat_props_.freundlich_n.data[i] = state.isotherm_freundlich_n_old[i][cell];
      alq_mat_props_.langmuir_b.data[i] = state.isotherm_langmuir_b_old[i][cell];
    }
  }

  // first order reaction rate cnst
  for (unsigned int i = 0; i != number_aqueous_kinetics_; ++i) {
    alq_mat_props_.aqueous_kinetic_rate_cnst.data[i] = state.aqueous_kinetics_rate[i][cell];
  }

  // Auxiliary data
  int num_aux_ints = chem_engine_->Sizes().num_aux_integers;
  int num_aux_doubles = chem_engine_->Sizes().num_aux_doubles;
  for (unsigned int i = 0; i != num_aux_ints; ++i) {
    alq_aux_data_.aux_ints.data[i] = (int)state.aux_data_old[i][cell];
  }
  for (unsigned int i = 0; i != num_aux_doubles; ++i) {
    alq_aux_data_.aux_doubles.data[i] = (int)state.aux_data_old[i+num_aux_ints][cell];
  }
}


/* *******************************************************************
*
******************************************************************* */
void
Alquimia_PK::CopyFromAlquimia(AlquimiaSubstate& state, const int cell, bool aux_output)
{
  for (unsigned i = 0; i != number_aqueous_components_; ++i) {
    state.aqueous_components[i][cell] = alq_state_.total_mobile.data[i];

    if (using_sorption_) {
      state.total_sorbed[i][cell] = alq_state_.total_immobile.data[i];
    }
  }

  // Free ion species.
  for (unsigned i = 0; i != number_aqueous_components_; ++i) {
    state.free_ion_species[i][cell] = aux_output.primary_free_ion_concentration.data[i];
  }

  // Mineral properties.
  for (unsigned i = 0; i != number_minerals_; ++i) {
    state.mineral_volume_fraction[i][cell] = alq_state_.mineral_volume_fraction.data[i];
    state.mineral_ssa[i][cell] = alq_state_.mineral_specific_surface_area.data[i];
    state.mineral_rate_constant[i][cell] = alq_mat_props_.mineral_rate_cnst.data[i];
  }

  // ion exchange
  for (unsigned unsigned i = 0; i != number_ion_exchange_sites_; ++i) {
    state.ion_exchange[i][cell] = alq_state_.cation_exchange_capacity.data[i];
  }

  // surface complexation
  for (unsigned unsigned i = 0; i != number_sorption_sites_; ++i) {
    state.sorption_sites[i][cell] = alq_state_.surface_site_density.data[i];
  }

  int num_aux_ints = chem_engine_->Sizes().num_aux_integers;
  int num_aux_doubles = chem_engine_->Sizes().num_aux_doubles;
  for (unsigned i = 0; i != num_aux_ints; ++i) {
    state.aux_data[i][cell] = (double)aux_data.aux_ints.data[i];
  }
  for (unsigned i = 0; i != num_aux_doubles; ++i) {
    state.aux_data[i+num_aux_its][cell] = aux_data.aux_doubles.data[i];
  }

  if (using_sorption_isotherms_) {
    for (unsigned unsigned i = 0; i != number_aqueous_components_; ++i) {
      state.isotherm_kd[i][cell] = alq_mat_props_.isotherm_kd.data[i];
      state.isotherm_freundlich_n[i][cell] = alq_mat_props_.freundlich_n.data[i];
      state.isotherm_langmuir_b[i][cell] = alq_mat_props_.langmuir_b.data[i];
    }
  }

  if (aux_output) {
    int numAqueousComplexes = chem_engine_->NumAqueousComplexes();

    
    for (int n = 0; n < map_[0].size(); ++n) (*aux_output_)[map_[0][n]][cell] = aux_output.pH;

    for (int n = 0; n < mineral_names_.size(); ++n)
      (*aux_output_)[map_[1][n]][cell] = aux_output.mineral_saturation_index.data[n];

    for (int n = 0; n < mineral_names_.size(); ++n)
      (*aux_output_)[map_[2][n]][cell] = aux_output.mineral_reaction_rate.data[n];

    for (int n = 0; n < primary_names_.size(); ++n)
      (*aux_output_)[map_[3][n]][cell] = aux_output.primary_free_ion_concentration.data[n];

    for (int n = 0; n < primary_names_.size(); ++n)
      (*aux_output_)[map_[4][n]][cell] = aux_output.primary_activity_coeff.data[n];

    for (int n = 0; n < numAqueousComplexes; ++n)
      (*aux_output_)[map_[5][n]][cell] = aux_output.secondary_free_ion_concentration.data[n];

    for (int n = 0; n < numAqueousComplexes; ++n)
      (*aux_output_)[map_[6][n]][cell] = aux_output.secondary_activity_coeff.data[n];
  }
}


/* *******************************************************************
* This helper advances the solution on a single cell within Amanzi's state.
* It returns the number of iterations taken to obtain the advanced solution,
* or -1 if an error occurred.
******************************************************************* */
int
Alquimia_PK::AdvanceSingleCell(double dt,
                               Teuchos::RCP<Epetra_MultiVector>& aqueous_components,
                               int cell)
{
  // Copy the state and property information from Amanzi's state within
  // this cell to Alquimia.
  //
  // NOTE: this should get set not to be hard-coded to Tags::DEFAULT, but
  // should use the same tag as transport.  See #673
  CopyToAlquimia(
    cell, aqueous_components, alq_mat_props_, alq_state_, alq_aux_data_, Tags::DEFAULT);

  int num_iterations = 0;
  if (alq_mat_props_.saturation > saturation_tolerance_) {
    bool success =
      chem_engine_->Advance(dt,
                            alq_mat_props_,
                            alq_state_,
                            alq_aux_data_,
                            alq_aux_output_,
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
  }

  // Move the information back into Amanzi's state, updating the given total concentration vector.
  CopyAlquimiaStateToAmanzi(
    cell, alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_, aqueous_components);

  return num_iterations;
}


/* *******************************************************************
* The MPC will call this function to signal to the process kernel that
* it has accepted the state update, thus, the PK should update
* possible auxilary state variables here
******************************************************************* */
void
Alquimia_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  Chemistry_PK::CommitStep(t_old, t_new, tag);

  // publish auxiliary data to state
  if (aux_output_ != Teuchos::null) {
    int counter = 0;
    for (int i = 0; i < aux_names_.size(); ++i) {
      auto& aux_state =
        *S_->GetW<CompositeVector>(aux_names_[i], tag_next_, passwd_).ViewComponent("cell");
      for (int j = 0; j < aux_subfield_names_[i].size(); ++j) {
        *aux_state(j) = *(*aux_output_)(counter++);
      }
    }
  }

}


/* *******************************************************************
* Auxiliary map from beacon to aux_names
******************************************************************* */
void
Alquimia_PK::InitializeAuxNamesMap_()
{
  int numAqueousComplexes = chem_engine_->NumAqueousComplexes();

  int counter = 0;
  map_.resize(7);
  for (int i = 0; i < aux_names_.size(); i++) {
    auto this_name = Keys::getVarName(aux_names_.at(i));
    if (this_name == "pH") {
      map_[0].push_back(counter++);
    } else if (this_name == "mineral_saturation_index") {
      // make sure all are present
      AMANZI_ASSERT(mineral_names_.size() == aux_subfield_names_[i].size());
      for (int j = 0; j < mineral_names_.size(); ++j) { map_[1].push_back(counter++); }
    } else if (this_name == "mineral_reaction_rate") {
      AMANZI_ASSERT(mineral_names_.size() == aux_subfield_names_[i].size());
      for (int j = 0; j < mineral_names_.size(); ++j) { map_[2].push_back(counter++); }
    } else if (this_name == "primary_free_ion_concentration") {
      AMANZI_ASSERT(primary_names_.size() == aux_subfield_names_[i].size());
      for (int j = 0; j < primary_names_.size(); ++j) { map_[3].push_back(counter++); }
    } else if (this_name == Keys::getVarName(primary_activity_coeff_key_)) {
      AMANZI_ASSERT(primary_names_.size() == aux_subfield_names_[i].size());
      for (int j = 0; j < primary_names_.size(); ++j) { map_[4].push_back(counter++); }
    } else if (this_name == "secondary_free_ion_concentration") {
      AMANZI_ASSERT(numAqueousComplexes == aux_subfield_names_[i].size());
      for (int j = 0; j < numAqueousComplexes; ++j) { map_[5].push_back(counter++); }
    } else if (this_name == Keys::getVarName(secondary_activity_coeff_key_)) {
      AMANZI_ASSERT(numAqueousComplexes == aux_subfield_names_[i].size());
      for (int j = 0; j < numAqueousComplexes; ++j) { map_[6].push_back(counter++); }
    }
  }
}

} // namespace AmanziChemistry
} // namespace Amanzi
