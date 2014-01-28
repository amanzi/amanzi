/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi Chemistry

License: see COPYRIGHT
Author: Jeffrey Johnson

This implements the Alquimia chemistry engine.

 ------------------------------------------------------------------------- */

#include <iostream>
#include <cstring>
#include <assert.h>
#include "ChemistryEngine.hh"
#include "errors.hh"
#include "exceptions.hh"

// Support for manipulating floating point exception handling.
#ifdef _GNU_SOURCE
#define AMANZI_USE_FENV
#include <fenv.h>
#endif

namespace Amanzi {
namespace AmanziChemistry {


ChemistryEngine::ChemistryEngine(const std::string& engineName, 
                                 const std::string& inputFile):
  chem_engine_name_(engineName),
  chem_engine_inputfile_(inputFile)
{
  Errors::Message msg;

  if (chem_engine_name_ != "PFloTran")
  {
    msg << "ChemistryEngine: Unsupported chemistry engine: '" << chem_engine_name_ << "'\n";
    msg << "  Options are 'PFlotran'.\n";
    Exceptions::amanzi_throw(msg);
  }

  // All alquimia function calls require a status object.
  AllocateAlquimiaEngineStatus(&chem_status_);
  CreateAlquimiaInterface(chem_engine_name_.c_str(), &chem_, &chem_status_);
  if (chem_status_.error != 0) 
  {
    std::cout << chem_status_.message << std::endl;
    msg << "ChemistryEngine: Could not create an interface to Alquimia.";
    Exceptions::amanzi_throw(msg); 
  }

  // Set up Alquimia, get sizes for data.
  chem_.Setup(chem_engine_inputfile_.c_str(),
              &engine_state_,
              &sizes_,
              &functionality_,
              &chem_status_);
  if (chem_status_.error != 0) 
  {
    std::cout << chem_status_.message << std::endl;
    PrintAlquimiaSizes(&sizes_);
    msg << "Error in creation of ChemistryEngine.";
    Exceptions::amanzi_throw(msg); 
  }

  // Allocate storage for additional Alquimia data.
  AllocateAlquimiaProblemMetaData(&sizes_, &chem_metadata_);

  // Get the problem metadata (species and mineral names, etc).
  chem_.GetProblemMetaData(&engine_state_,
                           &chem_metadata_,
                           &chem_status_);
  if (chem_status_.error != 0) 
  {
    std::cout << chem_status_.message << std::endl;
    PrintAlquimiaProblemMetaData(&chem_metadata_);
    msg << "Error in ChemistryEngine::Initialize";
    Exceptions::amanzi_throw(msg); 
  }
}

ChemistryEngine::~ChemistryEngine()
{
  chem_.Shutdown(&engine_state_, &chem_status_);
  FreeAlquimiaProblemMetaData(&chem_metadata_);

  // Delete the various geochemical conditions.
  for (std::map<std::string, AlquimiaGeochemicalCondition*>::iterator 
       iter = chem_conditions_.begin(); iter != chem_conditions_.end(); ++iter)
  {
    FreeAlquimiaGeochemicalCondition(iter->second);
    delete iter->second;
  }

  FreeAlquimiaEngineStatus(&chem_status_);
}

const std::string& ChemistryEngine::Name() const
{
  return chem_engine_name_;
}

bool ChemistryEngine::IsThreadSafe() const
{
  Errors::Message msg;

  if (not chem_initialized_)
  {
    msg << "ChemistryEngine: Cannot query before initialization!";
    Exceptions::amanzi_throw(msg); 
  }

  return functionality_.thread_safe;
}

int ChemistryEngine::NumPrimarySpecies() const
{
  return sizes_.num_primary;
}

int ChemistryEngine::NumAqueousComplexes() const
{
  return sizes_.num_aqueous_complexes;
}

int ChemistryEngine::NumSorbedSpecies() const
{
  return sizes_.num_sorbed;
}

void ChemistryEngine::GetPrimarySpeciesNames(std::vector<std::string>& speciesNames) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->primary_names.size;
  speciesNames.resize(N);
  for (int i = 0; i < N; ++i)
    speciesNames[i] = std::string(metadata->primary_names.data[i]);
}

int ChemistryEngine::NumMinerals() const
{
  return sizes_.num_kinetic_minerals;
}

void ChemistryEngine::GetMineralNames(std::vector<std::string>& mineralNames) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->mineral_names.size;
  mineralNames.resize(N);
  for (int i = 0; i < N; ++i)
    mineralNames[i] = std::string(metadata->mineral_names.data[i]);
}

int ChemistryEngine::NumSurfaceSites() const
{
  return sizes_.num_surface_sites;
}

void ChemistryEngine::GetSurfaceSiteNames(std::vector<std::string>& siteNames) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->surface_site_names.size;
  siteNames.resize(N);
  for (int i = 0; i < N; ++i)
    siteNames[i] = std::string(metadata->surface_site_names.data[i]);
}

int ChemistryEngine::NumIonExchangeSites() const
{
  return sizes_.num_ion_exchange_sites;
}

void ChemistryEngine::GetIonExchangeNames(std::vector<std::string>& ionExchangeNames) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->ion_exchange_names.size;
  ionExchangeNames.resize(N);
  for (int i = 0; i < N; ++i)
    ionExchangeNames[i] = std::string(metadata->ion_exchange_names.data[i]);
}

int ChemistryEngine::NumIsothermSpecies() const
{
  return sizes_.num_isotherm_species;
}

void ChemistryEngine::GetIsothermSpeciesNames(std::vector<std::string>& speciesNames) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->isotherm_species_names.size;
  speciesNames.resize(N);
  for (int i = 0; i < N; ++i)
    speciesNames[i] = std::string(metadata->isotherm_species_names.data[i]);
}

int ChemistryEngine::NumFreeIonSpecies() const
{
  return sizes_.num_primary;
}

void ChemistryEngine::EnforceCondition(const std::string& conditionName,
                                       const double time,
                                       const AlquimiaMaterialProperties& mat_props,
                                       AlquimiaState& chem_state,
                                       AlquimiaAuxiliaryData& aux_data,
                                       AlquimiaAuxiliaryOutputData& aux_output)
{
  Errors::Message msg;

  // Retrieve the chemical condition for the given name.
  std::map<std::string, AlquimiaGeochemicalCondition*>::iterator iter = chem_conditions_.find(conditionName);
  if (iter == chem_conditions_.end())
  {
    // NOTE: a condition with zero aqueous/mineral constraints is assumed to be defined in 
    // NOTE: the backend engine's input file. 
    AlquimiaGeochemicalCondition* condition = new AlquimiaGeochemicalCondition();
    int num_aq = 0, num_min = 0;
    AllocateAlquimiaGeochemicalCondition(kAlquimiaMaxStringLength, num_aq, num_min, condition);
    std::strcpy(condition->name, conditionName.c_str());

    // Add this to the conditions map.
    chem_conditions_[conditionName] = condition;
    iter = chem_conditions_.find(conditionName);
  }

#ifdef AMANZI_USE_FENV
  // Disable divide-by-zero floating point exceptions.
  int fpe_mask = fedisableexcept(FE_DIVBYZERO);
#endif 

  // Process the condition on the given array at the given time.
  // FIXME: Time is ignored for the moment.
  int ierr = 0;
  AlquimiaGeochemicalCondition* condition = iter->second;
  chem_.ProcessCondition(&engine_state_,
                         condition,
                         &(const_cast<AlquimiaMaterialProperties&>(mat_props)),
                         &chem_state,
                         &aux_data,
                         &chem_status_);

#ifdef AMANZI_USE_FENV
  // Re-enable pre-existing floating point exceptions.
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
    msg << "Error in enforcement of chemical condition '" << conditionName << "'";
    Exceptions::amanzi_throw(msg); 
  }  
#endif
}

void ChemistryEngine::Advance(const double delta_time,
                              const AlquimiaMaterialProperties& mat_props,
                              AlquimiaState& chem_state,
                              AlquimiaAuxiliaryData& aux_data,
                              AlquimiaAuxiliaryOutputData& aux_output,
                              int& num_iterations)
{
#ifdef AMANZI_USE_FENV
  // Disable divide-by-zero floating point exceptions.
  int fpe_mask = fedisableexcept(FE_DIVBYZERO);
#endif 

  // Advance the chemical reaction all operator-split-like.
  chem_.ReactionStepOperatorSplit(&engine_state_,
                                  (const_cast<double*>(&delta_time)),
                                  &(const_cast<AlquimiaMaterialProperties&>(mat_props)),
                                  &chem_state,
                                  &aux_data,
                                  &chem_status_);

#ifdef AMANZI_USE_FENV
  // Re-enable pre-existing floating point exceptions.
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
    msg << "Error in advance of chemical reactions.";
    Exceptions::amanzi_throw(msg); 
  }  
#endif

  // Write down the number of Newton iterations.
  num_iterations = chem_status_.num_newton_iterations;
}

const AlquimiaSizes& ChemistryEngine::Sizes() const
{
  return sizes_;
}

void ChemistryEngine::InitState(AlquimiaMaterialProperties& mat_props,
                                AlquimiaState& chem_state, 
                                AlquimiaAuxiliaryData& aux_data,
                                AlquimiaAuxiliaryOutputData& aux_output)
{
  AllocateAlquimiaMaterialProperties(&sizes_, &mat_props);
  AllocateAlquimiaState(&sizes_, &chem_state);
  AllocateAlquimiaAuxiliaryData(&sizes_, &aux_data);
  AllocateAlquimiaAuxiliaryOutputData(&sizes_, &aux_output);
}
                 
void ChemistryEngine::FreeState(AlquimiaMaterialProperties& mat_props,
                                AlquimiaState& chem_state,
                                AlquimiaAuxiliaryData& aux_data,
                                AlquimiaAuxiliaryOutputData& aux_output)
{
  FreeAlquimiaMaterialProperties(&mat_props);
  FreeAlquimiaState(&chem_state);
  FreeAlquimiaAuxiliaryData(&aux_data);
  FreeAlquimiaAuxiliaryOutputData(&aux_output);
}

} // namespace
} // namespace

