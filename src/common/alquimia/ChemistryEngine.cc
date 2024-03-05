/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Jeffrey Johnson
           Sergi Molins <smolins@lbl.gov>
*/

/*
  Alquimia

  This implements the Alquimia chemistry engine.
*/

#include <iostream>
#include <cstring>
#include <cstdio>
#include <assert.h>
#include "ChemistryEngine.hh"
#include "errors.hh"
#include "exceptions.hh"

// Support for manipulating floating point exception handling.
#ifdef _GNU_SOURCE
#  define AMANZI_USE_FENV
#  include <fenv.h>
#endif

namespace Amanzi {
namespace AmanziChemistry {

namespace {

void
CopyAlquimiaState(AlquimiaState* dest, AlquimiaState* src)
{
  dest->water_density = src->water_density;
  dest->porosity = src->porosity;
  dest->temperature = src->temperature;
  dest->aqueous_pressure = src->aqueous_pressure;
  memcpy(dest->total_mobile.data, src->total_mobile.data, sizeof(double) * src->total_mobile.size);
  memcpy(
    dest->total_immobile.data, src->total_immobile.data, sizeof(double) * src->total_immobile.size);
  memcpy(dest->mineral_volume_fraction.data,
         src->mineral_volume_fraction.data,
         sizeof(double) * src->mineral_volume_fraction.size);
  memcpy(dest->mineral_specific_surface_area.data,
         src->mineral_specific_surface_area.data,
         sizeof(double) * src->mineral_specific_surface_area.size);
  memcpy(dest->surface_site_density.data,
         src->surface_site_density.data,
         sizeof(double) * src->surface_site_density.size);
  memcpy(dest->cation_exchange_capacity.data,
         src->cation_exchange_capacity.data,
         sizeof(double) * src->cation_exchange_capacity.size);
}

// These functions are going into the next release of Alquimia.
void
CopyAlquimiaProperties(AlquimiaProperties* dest, AlquimiaProperties* src)
{
  dest->volume = src->saturation;
  dest->saturation = src->saturation;
  memcpy(dest->aqueous_kinetic_rate_cnst.data,
         src->aqueous_kinetic_rate_cnst.data,
         sizeof(double) * src->aqueous_kinetic_rate_cnst.size);
  memcpy(dest->mineral_rate_cnst.data,
         src->mineral_rate_cnst.data,
         sizeof(double) * src->mineral_rate_cnst.size);
  memcpy(dest->isotherm_kd.data, src->isotherm_kd.data, sizeof(double) * src->isotherm_kd.size);
  memcpy(dest->freundlich_n.data, src->freundlich_n.data, sizeof(double) * src->freundlich_n.size);
  memcpy(dest->langmuir_b.data, src->langmuir_b.data, sizeof(double) * src->langmuir_b.size);
}

void
CopyAlquimiaAuxiliaryData(AlquimiaAuxiliaryData* dest, AlquimiaAuxiliaryData* src)
{
  memcpy(dest->aux_ints.data, src->aux_ints.data, sizeof(int) * src->aux_ints.size);
  memcpy(dest->aux_doubles.data, src->aux_doubles.data, sizeof(double) * src->aux_doubles.size);
}

} // namespace


ChemistryEngine::ChemistryEngine(const std::string& engineName, const std::string& inputFile)
  : chem_engine_name_(engineName), chem_engine_inputfile_(inputFile)
{
  Errors::Message msg;

  // NOTE: Alquimia now has a "hands-off" mode in which Alquimia relies on
  // NOTE: reaction properties from the engine as opposed to the ones provided
  // NOTE: in Amanzi's input file. As stop-gap solution hands-off is indicated
  // NOTE: by the + sign at end of engine name
  bool hands_off = false;
  if (chem_engine_name_ == "PFloTran" || chem_engine_name_ == "CrunchFlow") hands_off = true;

  if (chem_engine_name_ != "PFloTran" && chem_engine_name_ != "CrunchFlow") {
    msg << "ChemistryEngine: Unsupported chemistry engine: '" << chem_engine_name_ << "'\n";
    msg << "  Options are 'PFloTran' or 'CrunchFlow'.\n";
    Exceptions::amanzi_throw(msg);
  }

  // All alquimia function calls require a status object.
  AllocateAlquimiaEngineStatus(&chem_status_);
  CreateAlquimiaInterface(chem_engine_name_.c_str(), &chem_, &chem_status_);
  if (chem_status_.error != 0) {
    std::cout << chem_status_.message << std::endl;
    msg << "ChemistryEngine: Could not create an interface to Alquimia.";
    Exceptions::amanzi_throw(msg);
  }

  // Set up Alquimia, get sizes for data.
  chem_.Setup(chem_engine_inputfile_.c_str(),
              hands_off,
              &engine_state_,
              &sizes_,
              &functionality_,
              &chem_status_);
  if (chem_status_.error != 0) {
    std::cout << chem_status_.message << std::endl;
    PrintAlquimiaSizes(&sizes_, stdout);
    msg << "Error in creation of ChemistryEngine.";
    Exceptions::amanzi_throw(msg);
  }

  // Allocate storage for additional Alquimia data.
  AllocateAlquimiaProblemMetaData(&sizes_, &chem_metadata_);

  // Get the problem metadata (species and mineral names, etc).
  chem_.GetProblemMetaData(&engine_state_, &chem_metadata_, &chem_status_);
  if (chem_status_.error != 0) {
    std::cout << chem_status_.message << std::endl;
    PrintAlquimiaProblemMetaData(&chem_metadata_, stdout);
    msg << "Error in ChemistryEngine::Initialize";
    Exceptions::amanzi_throw(msg);
  }
}

ChemistryEngine::~ChemistryEngine()
{
  chem_.Shutdown(&engine_state_, &chem_status_);
  FreeAlquimiaProblemMetaData(&chem_metadata_);

  // Delete the various geochemical conditions.
  for (GeochemicalConditionMap::iterator iter = chem_conditions_.begin();
       iter != chem_conditions_.end();
       ++iter) {
    FreeAlquimiaGeochemicalCondition(&iter->second->condition);
    FreeAlquimiaState(&iter->second->chem_state);
    FreeAlquimiaProperties(&iter->second->mat_props);
    FreeAlquimiaAuxiliaryData(&iter->second->aux_data);
    delete iter->second;
  }

  FreeAlquimiaEngineStatus(&chem_status_);
}

const std::string&
ChemistryEngine::Name() const
{
  return chem_engine_name_;
}

bool
ChemistryEngine::IsThreadSafe() const
{
  Errors::Message msg;

  if (not chem_initialized_) {
    msg << "ChemistryEngine: Cannot query before initialization!";
    Exceptions::amanzi_throw(msg);
  }

  return functionality_.thread_safe;
}

int
ChemistryEngine::NumPrimarySpecies() const
{
  return sizes_.num_primary;
}

int
ChemistryEngine::NumAqueousComplexes() const
{
  return sizes_.num_aqueous_complexes;
}

int
ChemistryEngine::NumSorbedSpecies() const
{
  return sizes_.num_sorbed;
}

void
ChemistryEngine::GetPrimarySpeciesNames(std::vector<std::string>& species_names) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->primary_names.size;
  species_names.resize(N);
  for (int i = 0; i < N; ++i) species_names[i] = std::string(metadata->primary_names.data[i]);
}

int
ChemistryEngine::NumMinerals() const
{
  return sizes_.num_minerals;
}

void
ChemistryEngine::GetMineralNames(std::vector<std::string>& mineral_names) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->mineral_names.size;
  mineral_names.resize(N);
  for (int i = 0; i < N; ++i) mineral_names[i] = std::string(metadata->mineral_names.data[i]);
}

int
ChemistryEngine::NumSurfaceSites() const
{
  return sizes_.num_surface_sites;
}

void
ChemistryEngine::GetSurfaceSiteNames(std::vector<std::string>& site_names) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->surface_site_names.size;
  site_names.resize(N);
  for (int i = 0; i < N; ++i) site_names[i] = std::string(metadata->surface_site_names.data[i]);
}

int
ChemistryEngine::NumIonExchangeSites() const
{
  return sizes_.num_ion_exchange_sites;
}

void
ChemistryEngine::GetIonExchangeNames(std::vector<std::string>& ion_exchange_names) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->ion_exchange_names.size;
  ion_exchange_names.resize(N);
  for (int i = 0; i < N; ++i)
    ion_exchange_names[i] = std::string(metadata->ion_exchange_names.data[i]);
}

int
ChemistryEngine::NumIsothermSpecies() const
{
  return sizes_.num_isotherm_species;
}

void
ChemistryEngine::GetIsothermSpeciesNames(std::vector<std::string>& species_names) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->isotherm_species_names.size;
  species_names.resize(N);
  for (int i = 0; i < N; ++i)
    species_names[i] = std::string(metadata->isotherm_species_names.data[i]);
}

int
ChemistryEngine::NumFreeIonSpecies() const
{
  return sizes_.num_primary;
}

void
ChemistryEngine::GetAuxiliaryOutputNames(
  std::vector<std::string>& aux_names,
  std::vector<std::vector<std::string>>& subfield_names) const
{
  aux_names.clear();
  subfield_names.clear();
  aux_names.emplace_back("pH");
  subfield_names.emplace_back(std::vector<std::string>{ "0" });

  const AlquimiaProblemMetaData* metadata = &chem_metadata_;

  // Mineral data -- one per mineral.
  int N = metadata->mineral_names.size;
  if (N > 0) {
    std::vector<std::string> mineral_names;
    for (int i = 0; i < N; ++i) { mineral_names.emplace_back(metadata->mineral_names.data[i]); }
    aux_names.emplace_back("mineral_saturation_index");
    subfield_names.emplace_back(mineral_names);
    aux_names.emplace_back("mineral_reaction_rate");
    subfield_names.emplace_back(std::move(mineral_names));
  }

  // Auxiliary data per primary species.
  N = metadata->primary_names.size;
  if (N > 0) {
    std::vector<std::string> primary_names;
    for (int i = 0; i < N; ++i) { primary_names.emplace_back(metadata->primary_names.data[i]); }
    aux_names.emplace_back("primary_free_ion_concentration");
    subfield_names.emplace_back(primary_names);
    aux_names.emplace_back("primary_activity_coeff");
    subfield_names.emplace_back(std::move(primary_names));
  }

  // Secondary auxiliary data.
  N = this->NumAqueousComplexes();
  if (N > 0) {
    std::vector<std::string> secondary_names;
    for (int i = 0; i < N; ++i) { secondary_names.emplace_back(std::to_string(i)); }
    aux_names.emplace_back("secondary_free_ion_concentration");
    subfield_names.emplace_back(secondary_names);
    aux_names.emplace_back("secondary_activity_coeff");
    subfield_names.emplace_back(std::move(secondary_names));
  }
}

int
ChemistryEngine::NumAqueousKinetics() const
{
  return sizes_.num_aqueous_kinetics;
}

void
ChemistryEngine::GetAqueousKineticNames(std::vector<std::string>& kinetics_names) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->aqueous_kinetic_names.size;
  kinetics_names.resize(N);
  for (int i = 0; i < N; ++i)
    kinetics_names[i] = std::string(metadata->aqueous_kinetic_names.data[i]);
}

void
ChemistryEngine::CreateCondition(const std::string& condition_name)
{
  // NOTE: a condition with zero aqueous/mineral constraints is assumed to be defined in
  // NOTE: the backend engine's input file.
  GeochemicalConditionData* condition = new GeochemicalConditionData();
  condition->processed = false;
  int num_aq = 0, num_min = 0;
  AllocateAlquimiaProperties(&sizes_, &condition->mat_props);
  AllocateAlquimiaGeochemicalCondition(
    kAlquimiaMaxStringLength, num_aq, num_min, &condition->condition);
  AllocateAlquimiaState(&sizes_, &condition->chem_state);
  AllocateAlquimiaAuxiliaryData(&sizes_, &condition->aux_data);
  std::strcpy(condition->condition.name, condition_name.c_str());

  // Add this to the conditions map.
  chem_conditions_[condition_name] = condition;
}


/* Mineral constraints will be discontinued in Alquimia -- see Sergi

void ChemistryEngine::AddMineralConstraint(const std::string& condition_name,
                                           const std::string& mineral_name,
                                           double volume_fraction,
                                           double specific_surface_area)
{
  assert(condition_name.length() > 0);
  assert(volume_fraction >= 0.0);
  assert(specific_surface_area >= 0.0);

  GeochemicalConditionMap::iterator iter = chem_conditions_.find(condition_name);
  if (iter != chem_conditions_.end())
  {
    AlquimiaGeochemicalCondition* condition = &iter->second->condition;

    // Do we have an existing constraint?
    int index = -1;
    for (int i = 0; i < condition->mineral_constraints.size; ++i)
    {
      if (!std::strcmp(condition->mineral_constraints.data[i].mineral_name, mineral_name.c_str()))
      {
        // Overwrite the old constraint.
        index = i;
        free(condition->mineral_constraints.data[index].mineral_name);
      }
    }
    if (index == -1)
    {
      // New constraint!
      index = condition->mineral_constraints.size;
      condition->mineral_constraints.size++;
      condition->mineral_constraints.data = (AlquimiaMineralConstraint*)realloc(condition->mineral_constraints.data, sizeof(AlquimiaMineralConstraint) * (index+1));
    }

    // Add the mineral constraint.
    condition->mineral_constraints.data[index].mineral_name = strdup(mineral_name.c_str());
    condition->mineral_constraints.data[index].volume_fraction = volume_fraction;
    condition->mineral_constraints.data[index].specific_surface_area = specific_surface_area;
  }
  else
  {
    Errors::Message msg;
    msg << "ChemistryEngine::AddMineralConstraint: no condition named '" << condition_name << "'.";
    Exceptions::amanzi_throw(msg);
  }
}
Mineral constraints will be discontinued in Alquimia -- see Sergi */


void
ChemistryEngine::AddAqueousConstraint(const std::string& condition_name,
                                      const std::string& primary_species_name,
                                      const std::string& constraint_type,
                                      const std::string& associated_species)
{
  assert(condition_name.length() > 0);
  assert(primary_species_name.length() > 0);
  assert((constraint_type == "total_aqueous") || (constraint_type == "charge") ||
         (constraint_type == "free") || (constraint_type == "mineral") ||
         (constraint_type == "gas") || (constraint_type == "pH"));

  GeochemicalConditionMap::iterator iter = chem_conditions_.find(condition_name);
  if (iter != chem_conditions_.end()) {
    AlquimiaGeochemicalCondition* condition = &iter->second->condition;

    /* Mineral constraints will be discontinued from Alquimia in the near future -- see Sergi

    // Is there a mineral constraint for the associated species?
    if (!associated_species.empty())
    {
      bool found_mineral = false;
      for (int i = 0; i < condition->mineral_constraints.size; ++i)
      {
        if (!std::strcmp(condition->mineral_constraints.data[i].mineral_name, associated_species.c_str()))
          found_mineral = true;
      }
      if (!found_mineral)
      {
        Errors::Message msg;
        msg << "ChemistryEngine::AddAqueousConstraint: the condition '" << condition_name << "' does not have a mineral constraint for '" << associated_species << "'.";
        Exceptions::amanzi_throw(msg);
      }
    }

    Mineral constraints will be discontinued from Alquimia in the near future -- see Sergi */

    // Do we have an existing constraint?
    int index = -1;
    for (int i = 0; i < condition->aqueous_constraints.size; ++i) {
      if (!std::strcmp(condition->aqueous_constraints.data[i].primary_species_name,
                       primary_species_name.c_str())) {
        // Overwrite the old constraint.
        index = i;
        free(condition->aqueous_constraints.data[index].primary_species_name);
        free(condition->aqueous_constraints.data[index].constraint_type);
        if (condition->aqueous_constraints.data[index].associated_species != NULL)
          free(condition->aqueous_constraints.data[index].associated_species);
      }
    }
    if (index == -1) {
      // New constraint!
      index = condition->aqueous_constraints.size;
      condition->aqueous_constraints.size++;
      condition->aqueous_constraints.data = (AlquimiaAqueousConstraint*)realloc(
        condition->aqueous_constraints.data, sizeof(AlquimiaAqueousConstraint) * (index + 1));
    }

    // Add the aqueous constraint.
    condition->aqueous_constraints.data[index].primary_species_name =
      strdup(primary_species_name.c_str());
    condition->aqueous_constraints.data[index].constraint_type = strdup(constraint_type.c_str());
    if (!associated_species.empty())
      condition->aqueous_constraints.data[index].associated_species =
        strdup(associated_species.c_str());
    else
      condition->aqueous_constraints.data[index].associated_species = NULL;
  } else {
    Errors::Message msg;
    msg << "ChemistryEngine::AddAqueousConstraint: no condition named '" << condition_name << "'.";
    Exceptions::amanzi_throw(msg);
  }
}

void
ChemistryEngine::EnforceCondition(const std::string& condition_name,
                                  const double time,
                                  const AlquimiaProperties& mat_props,
                                  AlquimiaState& chem_state,
                                  AlquimiaAuxiliaryData& aux_data,
                                  AlquimiaAuxiliaryOutputData& aux_output)
{
  Errors::Message msg;

  // Retrieve the chemical condition for the given name.
  GeochemicalConditionMap::iterator iter = chem_conditions_.find(condition_name);
  if (iter == chem_conditions_.end()) {
    CreateCondition(condition_name);
    iter = chem_conditions_.find(condition_name);
  }

#ifdef AMANZI_USE_FENV
  // Disable divide-by-zero floating point exceptions.
  int fpe_mask = fedisableexcept(FE_DIVBYZERO);
#endif

  AlquimiaGeochemicalCondition* condition = &iter->second->condition;
  AlquimiaProperties& nc_mat_props = const_cast<AlquimiaProperties&>(mat_props);
  if (!iter->second->processed) {
    // Copy the given state data into place for this condition.
    CopyAlquimiaProperties(&iter->second->mat_props, &nc_mat_props);
    CopyAlquimiaState(&iter->second->chem_state, &chem_state);
    CopyAlquimiaAuxiliaryData(&iter->second->aux_data, &aux_data);

    // Process the condition on the given array at the given time.
    // FIXME: Time is ignored for the moment.
    chem_.ProcessCondition(&engine_state_,
                           condition,
                           &iter->second->mat_props,
                           &iter->second->chem_state,
                           &iter->second->aux_data,
                           &chem_status_);
    iter->second->processed = true;
  }

  // Copy the constraint's data into place.
  // Here, condition names are used but Alquimia properties are
  // given as "Material" zones in Amanzi's inputthus disconnecting them
  // from "Initial" conditions (i.e. condition names)
  // For the sake of performance, we avoid here to re-process conditions
  // but we need to retain materical properties as provided in Amanzi inputs
  // CopyAlquimiaProperties(&nc_mat_props, &iter->second->mat_props);
  CopyAlquimiaState(&chem_state, &iter->second->chem_state);
  CopyAlquimiaAuxiliaryData(&aux_data, &iter->second->aux_data);

#ifdef AMANZI_USE_FENV
  // Re-enable pre-existing floating point exceptions.
  feclearexcept(fpe_mask);
  fpe_mask = feenableexcept(fpe_mask);
#endif

  // FIXME: Figure out a neutral parallel-friendly way to report errors.
  assert(chem_status_.error == 0);

#if 0
  if (chem_status_.error != 0)
    ierr = -1;

  // figure out if any of the processes threw an error, if so all processes will re-throw
  int recv = 0;
  mesh_->get_comm()->MaxAll(&ierr, &recv, 1);
  if (recv != 0)
  {
    msg << "Error in enforcement of chemical condition '" << condition_name << "'";
    Exceptions::amanzi_throw(msg);
  }
#endif
}

bool
ChemistryEngine::Advance(const double delta_time,
                         const AlquimiaProperties& mat_props,
                         AlquimiaState& chem_state,
                         AlquimiaAuxiliaryData& aux_data,
                         AlquimiaAuxiliaryOutputData& aux_output,
                         int& num_iterations, int natural_id)
{
#ifdef AMANZI_USE_FENV
  // Disable divide-by-zero floating point exceptions.
  int fpe_mask = fedisableexcept(FE_DIVBYZERO);
#endif

  // Advance the chemical reaction all operator-split-like.
  chem_.ReactionStepOperatorSplit(&engine_state_,
                                  delta_time,
                                  &(const_cast<AlquimiaProperties&>(mat_props)),
                                  &chem_state,
                                  &aux_data,
                                  natural_id,
                                  &chem_status_);

#ifdef AMANZI_USE_FENV
  // Re-enable pre-existing floating point exceptions.
  feclearexcept(fpe_mask);
  fpe_mask = feenableexcept(fpe_mask);
#endif

  // Retrieve auxiliary output.
  chem_.GetAuxiliaryOutput(&engine_state_,
                           &(const_cast<AlquimiaProperties&>(mat_props)),
                           &chem_state,
                           &aux_data,
                           &aux_output,
                           &chem_status_);

  // Did we succeed?
  if (chem_status_.error != kAlquimiaNoError) return false;

  // Did we converge?
  if (!chem_status_.converged) return false;

  // Write down the (maximum) number of Newton iterations.
  num_iterations = chem_status_.num_newton_iterations;
  return true;
}

const AlquimiaSizes&
ChemistryEngine::Sizes() const
{
  return sizes_;
}

void
ChemistryEngine::InitState(AlquimiaProperties& mat_props,
                           AlquimiaState& chem_state,
                           AlquimiaAuxiliaryData& aux_data,
                           AlquimiaAuxiliaryOutputData& aux_output)
{
  AllocateAlquimiaProperties(&sizes_, &mat_props);
  AllocateAlquimiaState(&sizes_, &chem_state);
  AllocateAlquimiaAuxiliaryData(&sizes_, &aux_data);
  AllocateAlquimiaAuxiliaryOutputData(&sizes_, &aux_output);

  // Make sure the auxiliary ints/doubles are zeroed out.
  std::fill(aux_data.aux_ints.data, aux_data.aux_ints.data + aux_data.aux_ints.size, 0);
  std::fill(aux_data.aux_doubles.data, aux_data.aux_doubles.data + aux_data.aux_doubles.size, 0.0);
}

void
ChemistryEngine::FreeState(AlquimiaProperties& mat_props,
                           AlquimiaState& chem_state,
                           AlquimiaAuxiliaryData& aux_data,
                           AlquimiaAuxiliaryOutputData& aux_output)
{
  FreeAlquimiaProperties(&mat_props);
  FreeAlquimiaState(&chem_state);
  FreeAlquimiaAuxiliaryData(&aux_data);
  FreeAlquimiaAuxiliaryOutputData(&aux_output);
}

} // namespace AmanziChemistry
} // namespace Amanzi
