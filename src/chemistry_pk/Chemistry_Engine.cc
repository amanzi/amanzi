/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi Chemistry

License: see COPYRIGHT
Author: Jeffrey Johnson

This implements the Alquimia chemistry engine.

 ------------------------------------------------------------------------- */

#include "Chemistry_Engine.hh"
#include "errors.hh"
#include "exceptions.hh"

// Support for manipulating floating point exception handling.
#ifdef _GNU_SOURCE
#define AMANZI_USE_FENV
#include <fenv.h>
#endif

namespace Amanzi {
namespace AmanziChemistry {


Chemistry_Engine::Chemistry_Engine(const Teuchos::ParameterList& param_list):
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

Chemistry_Engine::~Chemistry_Engine()
{
  if (chem_initialized_)
  {
    // Delete the context objects for the geochemical conditions.
    for (std::map<std::string, GeochemicalConditionContext*>::iterator 
         iter = chem_contexts_.begin(); iter != chem_contexts_.end(); ++iter)
    {
      delete iter->second;
    }

    chem_.Shutdown(&chem_data_.engine_state, &chem_status_);
    FreeAlquimiaData(&chem_data_);

    // Delete the various geochemical conditions.
    for (size_t i = 0; i < all_chem_conditions_.size(); ++i)
      FreeAlquimiaGeochemicalCondition(all_chem_conditions_[i]);

    FreeAlquimiaEngineStatus(&chem_status_);
  }
}

const std::string& Chemistry_Engine::Name() const
{
  return chem_engine_name_;
}

void Chemistry_Engine::Initialize()
{
  Errors::Message msg;

  // Read XML parameters from our input file.
  ReadXMLParameters();

  // We can't initialize chemistry twice in a row.
  if (chem_initialized_)
  {
    msg << "Alquimia_Chemistry_PK: initialized chemistry twice!";
    Exceptions::amanzi_throw(msg); 
  }

  // All alquimia function calls require a status object.
  AllocateAlquimiaEngineStatus(&chem_status_);
  CreateAlquimiaInterface(chem_engine_name_.c_str(), &chem_, &chem_status_);
  if (chem_status_.error != 0) 
  {
    std::cout << chem_status_.message << std::endl;
    msg << "Alquimia_Chemistry_PK: Could not create an interface to Alquimia.";
    Exceptions::amanzi_throw(msg); 
  }

  // Set up Alquimia, get sizes for data.
  chem_.Setup(chem_engine_inputfile_.c_str(),
              &chem_data_.engine_state,
              &chem_data_.sizes,
              &chem_data_.functionality,
              &chem_status_);
  if (chem_status_.error != 0) 
  {
    std::cout << chem_status_.message << std::endl;
    PrintAlquimiaSizes(&chem_data_.sizes);
    msg << "Error in Alquimia_Chemistry_PK::InitializeChemistry 1";
    Exceptions::amanzi_throw(msg); 
  }

  // Allocate data storage now that we have our space requirements.
  AllocateAlquimiaData(&chem_data_);

  // Get the problem metadata (species and mineral names, etc).
  chem_.GetProblemMetaData(&chem_data_.engine_state,
                           &chem_data_.meta_data,
                           &chem_status_);
  if (chem_status_.error != 0) 
  {
    std::cout << chem_status_.message << std::endl;
    PrintAlquimiaProblemMetaData(&chem_data_.meta_data);
    msg << "Error in Alquimia_Chemistry_PK::InitializeChemistry 1";
    Exceptions::amanzi_throw(msg); 
  }

  // FIXME: Make storage within Amanzi for auxiliary output.
  // FIXME: Use the Alquimia containers to do this.

  chem_initialized_ = true;
}

int Chemistry_Engine::NumPrimarySpecies() const
{
  return chem_data_.sizes.num_primary;
}

int Chemistry_Engine::NumAqueousComplexes() const
{
  return chem_data_.sizes.num_aqueous_complexes;
}

int Chemistry_Engine::NumSorbedSpecies() const
{
  return chem_data_.sizes.num_sorbed;
}

void Chemistry_Engine::GetPrimarySpeciesNames(std::vector<std::string>& speciesNames) const
{
  const AlquimiaProblemMetaData* metadata = &chem_data_.meta_data;
  int N = metadata->primary_names.size;
  speciesNames.resize(N);
  for (int i = 0; i < N; ++i)
    speciesNames[i] = std::string(metadata->primary_names.data[i]);
}

int Chemistry_Engine::NumMinerals() const
{
  return chem_data_.sizes.num_kinetic_minerals;
}

void Chemistry_Engine::GetMineralNames(std::vector<std::string>& mineralNames) const
{
  const AlquimiaProblemMetaData* metadata = &chem_data_.meta_data;
  int N = metadata->mineral_names.size;
  mineralNames.resize(N);
  for (int i = 0; i < N; ++i)
    mineralNames[i] = std::string(metadata->mineral_names.data[i]);
}

int Chemistry_Engine::NumSurfaceSites() const
{
  return chem_data_.sizes.num_surface_sites;
}

void Chemistry_Engine::GetSurfaceSiteNames(std::vector<std::string>& siteNames) const
{
  const AlquimiaProblemMetaData* metadata = &chem_data_.meta_data;
  int N = metadata->surface_site_names.size;
  siteNames.resize(N);
  for (int i = 0; i < N; ++i)
    siteNames[i] = std::string(metadata->surface_site_names.data[i]);
}

int Chemistry_Engine::NumIonExchangeSites() const
{
  return chem_data_.sizes.num_ion_exchange_sites;
}

void Chemistry_Engine::GetIonExchangeNames(std::vector<std::string>& ionExchangeNames) const
{
  const AlquimiaProblemMetaData* metadata = &chem_data_.meta_data;
  int N = metadata->ion_exchange_names.size;
  ionExchangeNames.resize(N);
  for (int i = 0; i < N; ++i)
    ionExchangeNames[i] = std::string(metadata->ion_exchange_names.data[i]);
}

int Chemistry_Engine::NumIsothermSpecies() const
{
  return chem_data_.sizes.num_isotherm_species;
}

void Chemistry_Engine::GetIsothermSpeciesNames(std::vector<std::string>& speciesNames) const
{
  const AlquimiaProblemMetaData* metadata = &chem_data_.meta_data;
  int N = metadata->isotherm_species_names.size;
  speciesNames.resize(N);
  for (int i = 0; i < N; ++i)
    speciesNames[i] = std::string(metadata->isotherm_species_names.data[i]);
}

int Chemistry_Engine::NumFreeIonSpecies() const
{
  return chem_data_.aux_output.primary_free_ion_concentration.size;
}

void Chemistry_Engine::EnforceCondition(const std::string& conditionName,
                                        double time,
                                        AlquimiaState* chem_state,
                                        AlquimiaMaterialProperties* mat_props,
                                        AlquimiaAuxiliaryData* aux_data,
                                        AlquimiaAuxiliaryOutputData* aux_output)
{
  Errors::Message msg;

  // Retrieve the chemical condition for the given name.
  std::map<std::string, AlquimiaGeochemicalCondition*>::iterator iter = chem_conditions_.find(conditionName);
  if (iter == chem_conditions_.end())
  {
    msg << "Chemistry_Engine::EnforceCondition: No condition '" << conditionName << "' was found!\n";
    Exceptions::amanzi_throw(msg); 
  }

  // Copy stuff into Alquimia.
  CopyIn(chem_state, mat_props, aux_data);

#ifdef AMANZI_USE_FENV
  // Disable divide-by-zero floating point exceptions.
  int fpe_mask = fedisableexcept(FE_DIVBYZERO);
#endif 

  // Process the condition on the given array at the given time.
  // FIXME: Time is ignored for the moment.
  int ierr = 0;
  AlquimiaGeochemicalCondition* condition = iter->second;
  chem_.ProcessCondition(&chem_data_.engine_state,
                         condition,
                         &chem_data_.material_properties,
                         &chem_data_.state,
                         &chem_data_.aux_data,
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

  // Copy stuff out again.
  CopyOut(chem_state, mat_props, aux_data, aux_output);
}

void Chemistry_Engine::Advance(double delta_time,
                               AlquimiaState* chem_state,
                               AlquimiaMaterialProperties* mat_props,
                               AlquimiaAuxiliaryData* aux_data,
                               AlquimiaAuxiliaryOutputData* aux_output,
                               int& num_iterations)
{
  CopyIn(chem_state, mat_props, aux_data);

#ifdef AMANZI_USE_FENV
  // Disable divide-by-zero floating point exceptions.
  int fpe_mask = fedisableexcept(FE_DIVBYZERO);
#endif 

  // Advance the chemical reaction all operator-split-like.
  chem_.ReactionStepOperatorSplit(&chem_data_.engine_state,
                                  &delta_time,
                                  &chem_data_.material_properties,
                                  &chem_data_.state,
                                  &chem_data_.aux_data,
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

  // Copy stuff out again.
  CopyOut(chem_state, mat_props, aux_data, aux_output);
}

GeochemicalConditionContext* Chemistry_Engine::ContextForCondition(const std::string& geochemical_condition_name)
{
  GeochemicalConditionContext* context = NULL;
  std::map<std::string, GeochemicalConditionContext*>::iterator iter = chem_contexts_.find(geochemical_condition_name);
  if (iter == chem_contexts_.end())
  {
    context = new GeochemicalConditionContext(this, geochemical_condition_name);
    chem_contexts_[geochemical_condition_name] = context;
  }
  else
    context = iter->second;
  return context;
}

void Chemistry_Engine::ParseChemicalConditions(const Teuchos::ParameterList& param_list,
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
      msg << "Chemistry_Engine::ParseChemicalConditions():\n";
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

void Chemistry_Engine::ReadXMLParameters() 
{
  Errors::Message msg;

  // Engine and engine input file.
  if (chem_param_list_.isParameter("Engine")) 
  {
    Errors::Message msg;
    chem_engine_name_ = chem_param_list_.get<std::string>("Engine");
    if (chem_engine_name_ != "PFloTran")
    {
      msg << "Chemistry_Engine::ReadXMLParameters(): \n";
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
//        message << "Chemistry_Engine::ReadXMLParameters(): unknown value in 'Auxiliary Data' list: " 
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

const AlquimiaSizes& Chemistry_Engine::Sizes() const
{
  return chem_data_.sizes;
}

AlquimiaState* Chemistry_Engine::NewState() const
{
  AlquimiaState* state = new AlquimiaState;
  AllocateAlquimiaState(&chem_data_.sizes, state);
  return state;
}

void Chemistry_Engine::DeleteState(AlquimiaState* state)
{
  FreeAlquimiaState(state);
  delete state;
}

AlquimiaMaterialProperties* Chemistry_Engine::NewMaterialProperties() const
{
  AlquimiaMaterialProperties* mat_props = new AlquimiaMaterialProperties;
  AllocateAlquimiaMaterialProperties(&chem_data_.sizes, mat_props);
  return mat_props;
}

void Chemistry_Engine::DeleteMaterialProperties(AlquimiaMaterialProperties* mat_props)
{
  FreeAlquimiaMaterialProperties(mat_props);
  delete mat_props;
}

AlquimiaAuxiliaryData* Chemistry_Engine::NewAuxiliaryData() const
{
  AlquimiaAuxiliaryData* aux_data = new AlquimiaAuxiliaryData;
  AllocateAlquimiaAuxiliaryData(&chem_data_.sizes, aux_data);
  return aux_data;
}

void Chemistry_Engine::DeleteAuxiliaryData(AlquimiaAuxiliaryData* aux_data)
{
  FreeAlquimiaAuxiliaryData(aux_data);
  delete aux_data;
}

AlquimiaAuxiliaryOutputData* Chemistry_Engine::NewAuxiliaryOutputData() const
{
  AlquimiaAuxiliaryOutputData* aux_output = new AlquimiaAuxiliaryOutputData;
  AllocateAlquimiaAuxiliaryOutputData(&chem_data_.sizes, aux_output);
  return aux_output;
}

void Chemistry_Engine::DeleteAuxiliaryOutputData(AlquimiaAuxiliaryOutputData* aux_output)
{
  FreeAlquimiaAuxiliaryOutputData(aux_output);
  delete aux_output;
}

void Chemistry_Engine::CopyIn(AlquimiaState* chem_state, AlquimiaMaterialProperties* mat_props, AlquimiaAuxiliaryData* aux_data)
{
  AlquimiaState* alquimia_state = &chem_data_.state;
  AlquimiaMaterialProperties* alquimia_mat_props = &chem_data_.material_properties;
  AlquimiaAuxiliaryData* alquimia_aux_data = &chem_data_.aux_data;

  // Chemical composition.
  for (int c = 0; c < alquimia_state->total_mobile.size; ++c)
    alquimia_state->total_mobile.data[c] = chem_state->total_mobile.data[c];
  for (int c = 0; c < alquimia_state->total_immobile.size; ++c)
    alquimia_state->total_immobile.data[c] = chem_state->total_immobile.data[c];

  // Mineral information.
  for (int m = 0; m < alquimia_state->mineral_volume_fraction.size; ++m) 
  {
    alquimia_state->mineral_volume_fraction.data[m] = chem_state->mineral_volume_fraction.data[m];
    alquimia_state->mineral_specific_surface_area.data[m] = chem_state->mineral_specific_surface_area.data[m];
  }

  // Ion exchange.
  for (int i = 0; i < alquimia_state->cation_exchange_capacity.size; ++i)
    alquimia_state->cation_exchange_capacity.data[i] = chem_state->cation_exchange_capacity.data[i];

  // Surface complexation.
  for (int s = 0; s < alquimia_state->surface_site_density.size; ++s)
    alquimia_state->surface_site_density.data[s] = chem_state->surface_site_density.data[s];

  // "Geometry."
  alquimia_state->water_density = chem_state->water_density; 
  alquimia_state->porosity = chem_state->porosity; 

  // Material properties.
  alquimia_mat_props->volume = mat_props->volume;
  alquimia_mat_props->saturation = mat_props->saturation;
  for (int i = 0; i < alquimia_mat_props->isotherm_kd.size; ++i)
    alquimia_mat_props->isotherm_kd.data[i] = mat_props->isotherm_kd.data[i];
  for (int i = 0; i < alquimia_mat_props->freundlich_n.size; ++i)
    alquimia_mat_props->freundlich_n.data[i] = mat_props->isotherm_kd.data[i];
  for (int i = 0; i < alquimia_mat_props->langmuir_b.size; ++i)
    alquimia_mat_props->langmuir_b.data[i] = mat_props->isotherm_kd.data[i];

  // Auxiliary data.
  for (int i = 0; i < alquimia_aux_data->aux_ints.size; ++i)
    alquimia_aux_data->aux_ints.data[i] = aux_data->aux_ints.data[i];
  for (int i = 0; i < alquimia_aux_data->aux_doubles.size; ++i)
    alquimia_aux_data->aux_doubles.data[i] = aux_data->aux_doubles.data[i];
}

void Chemistry_Engine::CopyOut(AlquimiaState* chem_state, AlquimiaMaterialProperties* mat_props, 
                               AlquimiaAuxiliaryData* aux_data, AlquimiaAuxiliaryOutputData* aux_output)
{
  AlquimiaState* alquimia_state = &chem_data_.state;
  AlquimiaMaterialProperties* alquimia_mat_props = &chem_data_.material_properties;
  AlquimiaAuxiliaryData* alquimia_aux_data = &chem_data_.aux_data;

  // Chemical composition.
  for (int c = 0; c < alquimia_state->total_mobile.size; ++c)
    chem_state->total_mobile.data[c] = alquimia_state->total_mobile.data[c];
  for (int c = 0; c < alquimia_state->total_immobile.size; ++c)
    chem_state->total_immobile.data[c] = alquimia_state->total_immobile.data[c];

  // Mineral information.
  for (int m = 0; m < alquimia_state->mineral_volume_fraction.size; ++m) 
  {
    chem_state->mineral_volume_fraction.data[m] = alquimia_state->mineral_volume_fraction.data[m];
    chem_state->mineral_specific_surface_area.data[m] = alquimia_state->mineral_specific_surface_area.data[m];
  }

  // Ion exchange.
  for (int i = 0; i < alquimia_state->cation_exchange_capacity.size; ++i)
    chem_state->cation_exchange_capacity.data[i] = alquimia_state->cation_exchange_capacity.data[i];

  // Surface complexation.
  for (int s = 0; s < alquimia_state->surface_site_density.size; ++s)
    chem_state->surface_site_density.data[s] = alquimia_state->surface_site_density.data[s];

  // "Geometry."
  chem_state->water_density = alquimia_state->water_density;
  chem_state->porosity = alquimia_state->porosity;

  // Material properties.
  mat_props->volume = alquimia_mat_props->volume;
  mat_props->saturation = alquimia_mat_props->saturation;
  for (int i = 0; i < alquimia_mat_props->isotherm_kd.size; ++i)
    mat_props->isotherm_kd.data[i] = alquimia_mat_props->isotherm_kd.data[i];
  for (int i = 0; i < alquimia_mat_props->freundlich_n.size; ++i)
    mat_props->isotherm_kd.data[i] = alquimia_mat_props->freundlich_n.data[i];
  for (int i = 0; i < alquimia_mat_props->langmuir_b.size; ++i)
    mat_props->isotherm_kd.data[i] = alquimia_mat_props->langmuir_b.data[i];

  // Auxiliary data.
  for (int i = 0; i < alquimia_aux_data->aux_ints.size; ++i)
    aux_data->aux_ints.data[i] = alquimia_aux_data->aux_ints.data[i];
  for (int i = 0; i < alquimia_aux_data->aux_doubles.size; ++i)
    aux_data->aux_doubles.data[i] = alquimia_aux_data->aux_doubles.data[i];

  // Auxiliary output data.
  if (aux_output != NULL)
  {
    const AlquimiaAuxiliaryOutputData* alquimia_aux_output = &chem_data_.aux_output;

    // Free ion concentrations.
    for (int i = 0; i < alquimia_aux_output->primary_free_ion_concentration.size; ++i)
      aux_output->primary_free_ion_concentration.data[i] = alquimia_aux_output->primary_free_ion_concentration.data[i];

    // Activity coefficients.
    for (int i = 0; i < alquimia_aux_output->primary_activity_coeff.size; ++i)
      aux_output->primary_activity_coeff.data[i] = alquimia_aux_output->primary_activity_coeff.data[i];
    for (int i = 0; i < alquimia_aux_output->secondary_activity_coeff.size; ++i)
      aux_output->secondary_activity_coeff.data[i] = alquimia_aux_output->secondary_activity_coeff.data[i];

    // pH.
    aux_output->pH = alquimia_aux_output->pH;
  }
}

// This subclass of Amanzi::Function provides an interface by which a geochemical condition can be 
// enforced on a given species.
class GeochemicalConcentrationFunction: public Function {

 public:

  // Constructs a GeochemicalConcentrationFunction that reports species concentrations 
  // computed by a GeochemicalConditionContext.
  explicit GeochemicalConcentrationFunction(int speciesIndex):
    index_(speciesIndex)
  {
  }

  // Destructor.
  ~GeochemicalConcentrationFunction();

  // Overridden methods.
  Function* clone() const
  {
    return new GeochemicalConcentrationFunction(index_);
  }

  double operator() (const double* xt) const
  {
    // No space/time dependence -- we just reach into our context for the answer.
    return value_;
  }

  double value_; // Updated by GeochemicalConditionContext.

 private:

  int index_;
};

GeochemicalConditionContext::GeochemicalConditionContext(Chemistry_Engine* chem_engine, 
                                                         const std::string& geochem_condition):
  chem_engine_(chem_engine_), condition_(geochem_condition), functions_(chem_engine->NumPrimarySpecies())
{
}

GeochemicalConditionContext::~GeochemicalConditionContext()
{
}

Teuchos::RCP<Function> GeochemicalConditionContext::speciesFunction(const std::string& species)
{
  // Find the index of the given species within the list kept in the chemistry engine.
  int speciesIndex = -1;
  std::vector<std::string> speciesNames;
  chem_engine_->GetPrimarySpeciesNames(speciesNames);
  for (size_t i = 0; i < speciesNames.size(); ++i)
  {
    if (species == speciesNames[i])
    {
      speciesIndex = (int)i;
      break;
    }
  }
  // FIXME: Check for speciesIndex == -1.
  GeochemicalConcentrationFunction* function = new GeochemicalConcentrationFunction(speciesIndex);
  functions_[speciesIndex] = Teuchos::RCP<Function>(function);
  return functions_[speciesIndex];
}

void GeochemicalConditionContext::EnforceCondition(double t,
                                                   AlquimiaState* chem_state,
                                                   AlquimiaMaterialProperties* mat_props,
                                                   AlquimiaAuxiliaryData* aux_data)
{
  // Enforce the condition.
  chem_engine_->EnforceCondition(condition_, t, chem_state, mat_props, aux_data);

  // Copy out the concentrations and broadcast them to our functions.
  int N = chem_engine_->NumPrimarySpecies();
  for (int i = 0; i < N; ++i)
  {
    Teuchos::RCP<GeochemicalConcentrationFunction> concFunc = functions_[i];
    concFunc->value_ = chem_state->total_mobile.data[i];
  }
}

} // namespace
} // namespace

