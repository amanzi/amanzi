/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi Chemistry

License: see COPYRIGHT
Author: Jeffrey Johnson

This implements the Alquimia chemistry engine.

 ------------------------------------------------------------------------- */

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


ChemistryEngine::ChemistryEngine(const Teuchos::ParameterList& param_list):
  main_param_list_(param_list), 
  chem_param_list_(),
  chem_initialized_(false),
  chem_engine_name_("PFloTran")
{
  assert(main_param_list_.isSublist("Chemistry"));
  chem_param_list_ = main_param_list_.sublist("Chemistry");

  // Set things up.
  Initialize();
}

ChemistryEngine::~ChemistryEngine()
{
  if (chem_initialized_)
  {
#if 0
    // Delete the context objects for the geochemical conditions.
    for (std::map<std::string, GeochemicalConditionContext*>::iterator 
         iter = chem_contexts_.begin(); iter != chem_contexts_.end(); ++iter)
    {
      delete iter->second;
    }
#endif

    chem_.Shutdown(engine_state_, &chem_status_);
    FreeAlquimiaProblemMetaData(&chem_metadata_);

    // Delete the various geochemical conditions.
    for (size_t i = 0; i < all_chem_conditions_.size(); ++i)
      FreeAlquimiaGeochemicalCondition(all_chem_conditions_[i]);

    FreeAlquimiaEngineStatus(&chem_status_);
  }
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

void ChemistryEngine::Initialize()
{
  Errors::Message msg;

  // Read XML parameters from our input file.
  ReadXMLParameters();

  // We can't initialize chemistry twice in a row.
  if (chem_initialized_)
  {
    msg << "ChemistryEngine: initialized chemistry twice!";
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
              engine_state_,
              &sizes_,
              &functionality_,
              &chem_status_);
  if (chem_status_.error != 0) 
  {
    std::cout << chem_status_.message << std::endl;
    PrintAlquimiaSizes(&sizes_);
    msg << "Error in ChemistryEngine::Initialize";
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

  chem_initialized_ = true;
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
                                       double time,
                                       AlquimiaState& chem_state,
                                       AlquimiaMaterialProperties& mat_props,
                                       AlquimiaAuxiliaryData& aux_data,
                                       AlquimiaAuxiliaryOutputData& aux_output)
{
  Errors::Message msg;

  // Retrieve the chemical condition for the given name.
  std::map<std::string, AlquimiaGeochemicalCondition*>::iterator iter = chem_conditions_.find(conditionName);
  if (iter == chem_conditions_.end())
  {
    msg << "ChemistryEngine::EnforceCondition: No condition '" << conditionName << "' was found!\n";
    Exceptions::amanzi_throw(msg); 
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
                         &mat_props,
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

void ChemistryEngine::Advance(double delta_time,
                              AlquimiaState& chem_state,
                              AlquimiaMaterialProperties& mat_props,
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
                                  &delta_time,
                                  &mat_props,
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

void ChemistryEngine::ParseChemicalConditions(const Teuchos::ParameterList& param_list,
                                              std::map<std::string, AlquimiaGeochemicalCondition*>& conditions)
{
  Errors::Message msg;

  // Go through the sublist containing the chemical conditions.
  for (Teuchos::ParameterList::ConstIterator iter = param_list.begin();
       iter != param_list.end(); ++iter)
  {
    // This parameter list contains sublists, each corresponding to a
    // condition. 
    std::string cond_name = param_list.name(iter);
    assert(param_list.isSublist(cond_name));
    const Teuchos::ParameterList& cond_sublist = param_list.sublist(cond_name);
 
    // At the moment, we only support constraints supplied through the backend engine.
    std::string engine_constraint = chem_engine_name_ + std::string(" Constraint");
    if (!cond_sublist.isType<std::string>(engine_constraint))
    {
      msg << "ChemistryEngine::ParseChemicalConditions():\n";
      msg << "  Geochemical condition '" << cond_name << "' has no '" << engine_constraint << "' entry.\n";
      msg << "  (Must be of type Array(string).)\n";
      Exceptions::amanzi_throw(msg);
    }
    std::string constraint_name = cond_sublist.get<std::string>(engine_constraint);

    // NOTE: a condition with zero aqueous/mineral constraints is assumed to be defined in 
    // NOTE: the backend engine's input file. 
    AlquimiaGeochemicalCondition* condition = new AlquimiaGeochemicalCondition();
    int num_aq = 0, num_min = 0;
    AllocateAlquimiaGeochemicalCondition(kAlquimiaMaxStringLength, num_aq, num_min, condition);
    strcpy(condition->name, constraint_name.c_str());

    // Add this to the conditions map.
    conditions[cond_name] = condition;
  }
}

void ChemistryEngine::ReadXMLParameters() 
{
  Errors::Message msg;

  // Engine and engine input file.
  if (chem_param_list_.isParameter("Engine")) 
  {
    Errors::Message msg;
    chem_engine_name_ = chem_param_list_.get<std::string>("Engine");
    if (chem_engine_name_ != "PFloTran")
    {
      msg << "ChemistryEngine::ReadXMLParameters(): \n";
      msg << "  Unsupported chemistry engine: '" << chem_engine_name_ << "'\n";
      msg << "  Options are 'PFlotran'.\n";
      Exceptions::amanzi_throw(msg);
    }
  }
  if (chem_param_list_.isParameter("Engine Input File"))
  {
    chem_engine_inputfile_ = chem_param_list_.get<std::string>("Engine Input File");
    // FIXME: Should we try to access the file here?
  }

  // --------------------------------------------------------------------------
  //
  // Auxiliary Data
  //
  // --------------------------------------------------------------------------
  aux_names_.clear();
  if (chem_param_list_.isParameter("Auxiliary Data")) 
  {
    Teuchos::Array<std::string> names = chem_param_list_.get<Teuchos::Array<std::string> >("Auxiliary Data");
    for (Teuchos::Array<std::string>::const_iterator name = names.begin();
         name != names.end(); ++name) 
    {
      if (*name == "pH") 
      {
        aux_names_.push_back(*name);
      } 
//      else 
//      {
//        std::stringstream message;
//        message << "ChemistryEngine::ReadXMLParameters(): unknown value in 'Auxiliary Data' list: " 
//                << *name << std::endl;
//        chem_out->Write(kWarning, message);
//      }
    }
  }

  // Now set up chemical conditions based on initial and boundary conditions
  // elsewhere in the file.
  assert(chem_conditions_.empty());
  assert(all_chem_conditions_.empty());

  // Initial conditions.
  if (!main_param_list_.isSublist("state"))
//  if (!main_param_list_.isSublist("Initial Conditions"))
  {
    msg << "Chemistry_PK::ReadXMLParameters(): \n";
    msg << "  No 'State' sublist was found!\n";
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::ParameterList state_list = main_param_list_.sublist("state");
  Teuchos::ParameterList initial_conditions = state_list.sublist("initial conditions");
  std::map<std::string, AlquimiaGeochemicalCondition*> init_cond;
  ParseChemicalConditions(initial_conditions, init_cond);
  if (init_cond.empty())
  {
    msg << "Chemistry_PK::ReadXMLParameters(): \n";
    msg << "  No chemical conditions were found in 'initial conditions'!\n";
    Exceptions::amanzi_throw(msg);
  }

  // Boundary conditions.
  if (!chem_param_list_.isSublist("Boundary Conditions"))
  {
    msg << "Chemistry_PK::ReadXMLParameters(): \n";
    msg << "  No boundary conditions were found!\n";
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::ParameterList boundary_conditions = chem_param_list_.sublist("Boundary Conditions");
  std::map<std::string, AlquimiaGeochemicalCondition*> bound_cond;
  ParseChemicalConditions(boundary_conditions, bound_cond);
  if (bound_cond.empty())
  {
    msg << "Chemistry_PK::ReadXMLParameters(): \n";
    msg << "  No chemical boundary conditions were found!\n";
    Exceptions::amanzi_throw(msg);
  }

  // Aggregate the initial and boundary conditions into one big table.
  for (std::map<std::string, AlquimiaGeochemicalCondition*>::const_iterator 
       iter = init_cond.begin(); iter != init_cond.end(); ++iter)
  {
    chem_conditions_[iter->first] = iter->second;
    all_chem_conditions_.push_back(iter->second);
  }
  for (std::map<std::string, AlquimiaGeochemicalCondition*>::const_iterator 
       iter = bound_cond.begin(); iter != bound_cond.end(); ++iter)
  {
    // Make sure there's no initial condition with the same name as this boundary condition.
    if (chem_conditions_.find(iter->first) != chem_conditions_.end())
    {
      msg << "Chemistry_PK::ReadXMLParameters(): \n";
      msg << "  Chemical boundary condition '" << iter->first << "' has the same name as an initial condition!\n";
      Exceptions::amanzi_throw(msg);
    }
    chem_conditions_[iter->first] = iter->second;
    all_chem_conditions_.push_back(iter->second);
  }
  // FIXME: Need to revisit the purpose of the ownership of chemical conditions in 
  // FIXME: all_chem_conditions_.
       
}  // end ReadXMLParameters()

const AlquimiaSizes& ChemistryEngine::Sizes() const
{
  return sizes_;
}

void ChemistryEngine::InitState(AlquimiaState& chem_state, 
                                AlquimiaMaterialProperties& mat_props,
                                AlquimiaAuxiliaryData& aux_data,
                                AlquimiaAuxiliaryOutputData& aux_output)
{
  AllocateAlquimiaState(&sizes_, &chem_state);
  AllocateAlquimiaMaterialProperties(&sizes_, &mat_props);
  AllocateAlquimiaAuxiliaryData(&sizes_, &aux_data);
  AllocateAlquimiaAuxiliaryOutputData(&sizes_, &aux_output);
  FreeAlquimiaAuxiliaryOutputData(&aux_output);
}
                 
void ChemistryEngine::FreeState(AlquimiaState& chem_state,
                                AlquimiaMaterialProperties& mat_props,
                                AlquimiaAuxiliaryData& aux_data,
                                AlquimiaAuxiliaryOutputData& aux_output)
{
  FreeAlquimiaState(&chem_state);
  FreeAlquimiaMaterialProperties(&mat_props);
  FreeAlquimiaAuxiliaryData(&aux_data);
}

} // namespace
} // namespace

