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

void Chemistry_Engine::CopyFromBuffersToAlquimia(const double* component_concentrations,
                                                 const double* sorbed_components,
                                                 const double* mineral_volume_fractions,
                                                 const double* mineral_specific_surface_areas,
                                                 const double* cation_exchange_capacity,
                                                 const double* sorption_sites,
                                                 double water_density,
                                                 double porosity,
                                                 double volume,
                                                 double saturation,
                                                 const double* isotherm_kd,
                                                 const double* isotherm_freundlich_n,
                                                 const double* isotherm_langmuir_b)
{
  AlquimiaState* alquimia_state = &chem_data_.state;
  AlquimiaMaterialProperties* alquimia_mat_props = &chem_data_.material_properties;
  
  // Chemical composition.
  for (int c = 0; c < alquimia_state->total_mobile.size; ++c)
  {
    alquimia_state->total_mobile.data[c] = component_concentrations[c];
    if (this->NumSorbedSpecies() > 0)
      alquimia_state->total_immobile.data[c] = sorbed_components[c];
  }

  // Mineral information.
  for (int m = 0; m < alquimia_state->mineral_volume_fraction.size; ++m) 
  {
    alquimia_state->mineral_volume_fraction.data[m] = mineral_volume_fractions[m];
    alquimia_state->mineral_specific_surface_area.data[m] = mineral_specific_surface_areas[m];
  }

  // Ion exchange.
  for (int i = 0; i < alquimia_state->cation_exchange_capacity.size; ++i)
    alquimia_state->cation_exchange_capacity.data[i] = cation_exchange_capacity[i];

  // Surface complexation.
  for (int s = 0; s < alquimia_state->surface_site_density.size; ++s)
    alquimia_state->surface_site_density.data[s] = sorption_sites[s];

  // "Geometry."
  alquimia_state->water_density = water_density; 
  alquimia_state->porosity = porosity; 

  // Material properties.
  alquimia_mat_props->volume = volume;
  alquimia_mat_props->saturation = saturation;
  for (int i = 0; i < alquimia_mat_props->isotherm_kd.size; ++i)
    alquimia_mat_props->isotherm_kd.data[i] = isotherm_kd[i];
  for (int i = 0; i < alquimia_mat_props->freundlich_n.size; ++i)
    alquimia_mat_props->freundlich_n.data[i] = isotherm_kd[i];
  for (int i = 0; i < alquimia_mat_props->langmuir_b.size; ++i)
    alquimia_mat_props->langmuir_b.data[i] = isotherm_kd[i];
}

void Chemistry_Engine::CopyFromAlquimiaToBuffers(double* component_concentrations,
                                                 double* sorbed_components,
                                                 double* mineral_volume_fractions,
                                                 double* mineral_specific_surface_areas,
                                                 double* cation_exchange_capacity,
                                                 double* sorption_sites,
                                                 double& water_density,
                                                 double& porosity,
                                                 double& volume,
                                                 double& saturation,
                                                 double* isotherm_kd,
                                                 double* isotherm_freundlich_n,
                                                 double* isotherm_langmuir_b,
                                                 double* free_ion_species,
                                                 double* primary_activity_coeffs,
                                                 double* secondary_activity_coeffs,
                                                 double* ion_exchange_ref_cation_concs,
                                                 double* surface_complex_free_site_concs,
                                                 double& pH) const
{
  const AlquimiaState* alquimia_state = &chem_data_.state;

  // Chemical composition.
  for (int c = 0; c < alquimia_state->total_mobile.size; ++c)
  {
    component_concentrations[c] = alquimia_state->total_mobile.data[c];
    if (this->NumSorbedSpecies() > 0)
      sorbed_components[c] = alquimia_state->total_immobile.data[c];
  }

  // Mineral information.
  for (int m = 0; m < alquimia_state->mineral_volume_fraction.size; ++m) 
  {
    mineral_volume_fractions[m] = alquimia_state->mineral_volume_fraction.data[m];
    mineral_specific_surface_areas[m] = alquimia_state->mineral_specific_surface_area.data[m];
  }

  // Ion exchange.
  for (int i = 0; i < alquimia_state->cation_exchange_capacity.size; ++i)
    cation_exchange_capacity[i] = alquimia_state->cation_exchange_capacity.data[i];

  // Surface complexation.
  for (int s = 0; s < alquimia_state->surface_site_density.size; ++s)
    sorption_sites[s] = alquimia_state->surface_site_density.data[s];

  // "Geometry."
  water_density = alquimia_state->water_density;
  porosity = alquimia_state->porosity;

  // Material properties.
  const AlquimiaMaterialProperties* alquimia_mat_props = &chem_data_.material_properties;
  volume = alquimia_mat_props->volume;
  saturation = alquimia_mat_props->saturation;
  for (int i = 0; i < alquimia_mat_props->isotherm_kd.size; ++i)
    isotherm_kd[i] = alquimia_mat_props->isotherm_kd.data[i];
  for (int i = 0; i < alquimia_mat_props->freundlich_n.size; ++i)
    isotherm_kd[i] = alquimia_mat_props->freundlich_n.data[i];
  for (int i = 0; i < alquimia_mat_props->langmuir_b.size; ++i)
    isotherm_kd[i] = alquimia_mat_props->langmuir_b.data[i];
  
  // Auxiliary data follows.
  const AlquimiaAuxiliaryOutputData* aux_output = &chem_data_.aux_output;

  // Free ion concentrations.
  for (int i = 0; i < aux_output->primary_free_ion_concentration.size; ++i)
    free_ion_species[i] = aux_output->primary_free_ion_concentration.data[i];

  // Activity coefficients.
  for (int i = 0; i < aux_output->primary_activity_coeff.size; ++i)
    primary_activity_coeffs[i] = aux_output->primary_activity_coeff.data[i];
  for (int i = 0; i < aux_output->secondary_activity_coeff.size; ++i)
    secondary_activity_coeffs[i] = aux_output->secondary_activity_coeff.data[i];

  // pH.
  pH = aux_output->pH;

  // The ion exchange reference cation sorbed concentration and surface complexation 
  // free site concentrations are not made available in an engine-neutral format and 
  // are not implemented currently. We can do this for the PFlotran engine, but 
  // this doesn't seem like a good approach looking forward.
#if 0
  // The layout of auxiliary 
  //   free ion conc <N_primary>
  //   primary_activity_coeff <N_primary>
  //   secondary_activity_coeff <N_aqueous_complexes>
  //   ion exchange ref cation sorbed conc <N_ion_exchange_sites>
  //   surface complexation free site conc <N_surface_sites>

  for (int i = 0; i < this->NumPrimarySpecies(); ++i, ++offset) 
    free_ion_species[i] = chem_data_.aux_data.aux_doubles.data[offset];

  // Ion exchange ref cation concentrations.
  for (int i = 0; i < this->NumIonExchangeSites(); ++i, ++offset) 
    ion_exchange_ref_cation_concs[i] = chem_data_.aux_data.aux_doubles.data[offset];

  // Surface complexation.
  if (this->NumSorbedSpecies() > 0)
  {
    for (int i = 0; i < this->NumSurfaceSites(); ++i, ++offset) 
      surface_complex_free_site_concs[i] = chem_data_.aux_data.aux_doubles.data[offset];
  }
#endif
}

void Chemistry_Engine::EnforceCondition(const std::string& conditionName,
                                        double time,
                                        double* component_concentrations,
                                        double* sorbed_components,
                                        double* mineral_volume_fractions,
                                        double* mineral_specific_surface_areas,
                                        double* cation_exchange_capacity,
                                        double* sorption_sites,
                                        double& water_density,
                                        double& porosity,
                                        double& volume,
                                        double& saturation,
                                        double* isotherm_kd,
                                        double* isotherm_freundlich_n,
                                        double* isotherm_langmuir_b,
                                        double* free_ion_species,
                                        double* primary_activity_coeffs,
                                        double* secondary_activity_coeffs,
                                        double* ion_exchange_ref_cation_concs,
                                        double* surface_complex_free_site_concs,
                                        double& pH)
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
  this->CopyFromBuffersToAlquimia(component_concentrations,
                                  sorbed_components,
                                  mineral_volume_fractions,
                                  mineral_specific_surface_areas,
                                  cation_exchange_capacity,
                                  sorption_sites,
                                  water_density,
                                  porosity,
                                  volume,
                                  saturation,
                                  isotherm_kd,
                                  isotherm_freundlich_n,
                                  isotherm_langmuir_b);

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
  this->CopyFromAlquimiaToBuffers(component_concentrations,
                                  sorbed_components,
                                  mineral_volume_fractions,
                                  mineral_specific_surface_areas,
                                  cation_exchange_capacity,
                                  sorption_sites,
                                  water_density,
                                  porosity,
                                  volume,
                                  saturation,
                                  isotherm_kd,
                                  isotherm_freundlich_n,
                                  isotherm_langmuir_b,
                                  free_ion_species,
                                  primary_activity_coeffs,
                                  secondary_activity_coeffs,
                                  ion_exchange_ref_cation_concs,
                                  surface_complex_free_site_concs,
                                  pH);
}

void Chemistry_Engine::Advance(double delta_time,
                               double* component_concentrations,
                               double* sorbed_components,
                               double* mineral_volume_fractions,
                               double* mineral_specific_surface_areas,
                               double* cation_exchange_capacity,
                               double* sorption_sites,
                               double& water_density,
                               double& porosity,
                               double& volume,
                               double& saturation,
                               double* isotherm_kd,
                               double* isotherm_freundlich_n,
                               double* isotherm_langmuir_b,
                               double* free_ion_species,
                               double* primary_activity_coeffs,
                               double* secondary_activity_coeffs,
                               double* ion_exchange_ref_cation_concs,
                               double* surface_complex_free_site_concs,
                               double& pH,
                               int& num_iterations)
{
  // Copy stuff into Alquimia.
  this->CopyFromBuffersToAlquimia(component_concentrations,
                                  sorbed_components,
                                  mineral_volume_fractions,
                                  mineral_specific_surface_areas,
                                  cation_exchange_capacity,
                                  sorption_sites,
                                  water_density,
                                  porosity,
                                  volume,
                                  saturation,
                                  isotherm_kd,
                                  isotherm_freundlich_n,
                                  isotherm_langmuir_b);

  // Advance the chemical reaction all operator-split-like.
  chem_.ReactionStepOperatorSplit(&chem_data_.engine_state,
                                  &delta_time,
                                  &chem_data_.material_properties,
                                  &chem_data_.state,
                                  &chem_data_.aux_data,
                                  &chem_status_);

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
  this->CopyFromAlquimiaToBuffers(component_concentrations,
                                  sorbed_components,
                                  mineral_volume_fractions,
                                  mineral_specific_surface_areas,
                                  cation_exchange_capacity,
                                  sorption_sites,
                                  water_density,
                                  porosity,
                                  volume,
                                  saturation,
                                  isotherm_kd,
                                  isotherm_freundlich_n,
                                  isotherm_langmuir_b,
                                  free_ion_species,
                                  primary_activity_coeffs,
                                  secondary_activity_coeffs,
                                  ion_exchange_ref_cation_concs,
                                  surface_complex_free_site_concs,
                                  pH);
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


} // namespace
} // namespace

