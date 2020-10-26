/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

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
                         const Teuchos::RCP<TreeVector>& soln) :
  Chemistry_PK(),
  soln_(soln),
  max_time_step_(9.9e9),
  chem_initialized_(false),
  current_time_(0.0),
  saved_time_(0.0)
{
  S_ = S;
  glist_ = glist;

  // extract pk name
  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // create pointer to the chemistry parameter list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  cp_list_ = Teuchos::sublist(pk_list, pk_name, true);

  domain_ = cp_list_->get<std::string>("domain name", "domain");

  // obtain key of fields
  tcc_key_ = Keys::readKey(*cp_list_,domain_, "total component concentration", "total_component_concentration"); 

  poro_key_ = Keys::readKey(*cp_list_, domain_, "porosity", "porosity");
  saturation_key_ = Keys::readKey(*cp_list_, domain_, "saturation liquid", "saturation_liquid");
  fluid_den_key_ = Keys::readKey(*cp_list_, domain_, "mass density liquid", "mass_density_liquid");
  
  min_vol_frac_key_ = Keys::readKey(*cp_list_, domain_, "mineral volume fractions", "mineral_volume_fractions");
  min_ssa_key_ = Keys::readKey(*cp_list_, domain_, "mineral specific surface area", "mineral_specific_surface_area");
  sorp_sites_key_ = Keys::readKey(*cp_list_, domain_, "sorption sites", "sorption_sites");
  surf_cfsc_key_ = Keys::readKey(*cp_list_, domain_, "surface complex free site conc", "surface_complex_free_site_conc");
  total_sorbed_key_ = Keys::readKey(*cp_list_, domain_, "total sorbed", "total_sorbed");
  isotherm_kd_key_ = Keys::readKey(*cp_list_, domain_, "isotherm_kd", "isotherm_kd");
  isotherm_freundlich_n_key_ = Keys::readKey(*cp_list_, domain_, "isotherm freundlich_n", "isotherm_freundlich_n");
  isotherm_langmuir_b_key_ = Keys::readKey(*cp_list_, domain_, "isotherm langmuir_b", "isotherm_langmuir_b");
  free_ion_species_key_ = Keys::readKey(*cp_list_, domain_, "free ion species", "free_ion_species");
  primary_activity_coeff_key_ = Keys::readKey(*cp_list_, domain_, "primary activity coeff", "primary_activity_coeff");

  ion_exchange_sites_key_ = Keys::readKey(*cp_list_, domain_, "ion exchange sites", "ion_exchange_sites");
  ion_exchange_ref_cation_conc_key_ = Keys::readKey(*cp_list_, domain_, "ion exchange ref cation conc", "ion_exchange_ref_cation_conc");
  secondary_activity_coeff_key_ = Keys::readKey(*cp_list_, domain_, "secondary activity coeff", "secondary_activity_coeff");
  alquimia_aux_data_key_ = Keys::readKey(*cp_list_, domain_, "alquimia aux data", "alquimia_aux_data");

  ion_exchange_ref_cation_conc_key_ = Keys::readKey(*cp_list_, domain_,"ion exchange ref cation conc", "ion_exchange_ref_cation_conc");
  secondary_activity_coeff_key_ = Keys::readKey(*cp_list_, domain_,"secondary activity coeff", "secondary_activity_coeff");
  alquimia_aux_data_key_ = Keys::readKey(*cp_list_, domain_,"alquimia aux data", "alquimia_aux_data");
  mineral_rate_constant_key_ = Keys::readKey(*cp_list_, domain_,"mineral rate constant", "mineral_rate_constant");
  first_order_decay_constant_key_ = Keys::readKey(*cp_list_, domain_,"first order decay constant", "first_order_decay_constant");


  // collect high-level information about the problem
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(glist, "state", true);

  InitializeMinerals(cp_list_);
  if (cp_list_->isSublist("initial conditions")) {
    // ATS-style input spec -- initial conditions in the PK
    InitializeSorptionSites(cp_list_, cp_list_);
  } else {
    // Amanzi-style input spec -- initial conditions in State
    InitializeSorptionSites(cp_list_, state_list);
  }

  // create chemistry engine. (should we do it later in Setup()?)
  if (!cp_list_->isParameter("engine")) {
    Errors::Message msg;
    msg << "No 'engine' parameter found in the parameter list for 'Chemistry'.\n";
    Exceptions::amanzi_throw(msg);
  }
  if (!cp_list_->isParameter("engine input file")) {
    Errors::Message msg;
    msg << "No 'engine input file' parameter found in the parameter list for 'Chemistry'.\n";
    Exceptions::amanzi_throw(msg);
  }
  std::string engine_name = cp_list_->get<std::string>("engine");
  std::string engine_inputfile = cp_list_->get<std::string>("engine input file");
  chem_engine_ = Teuchos::rcp(new AmanziChemistry::ChemistryEngine(engine_name, engine_inputfile));

  // grab the component names
  comp_names_.clear();
  chem_engine_->GetPrimarySpeciesNames(comp_names_);

  number_aqueous_components_ = comp_names_.size();
  number_free_ion_ = number_aqueous_components_;
  number_total_sorbed_ = number_aqueous_components_;

  chem_engine_->GetAqueousKineticNames(aqueous_kinetics_names_);
  number_aqueous_kinetics_ = aqueous_kinetics_names_.size();
  
  // verbosity object
  vo_ = Teuchos::rcp(new VerboseObject("Alquimia_PK:" + domain_, *cp_list_));
}


/* *******************************************************************
* Destroy ansilary data structures.
******************************************************************* */
Alquimia_PK::~Alquimia_PK() 
{
  if (chem_initialized_)
    chem_engine_->FreeState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);
}


/* ******************************************************************
* Register fields and evaluators with the State
******************************************************************* */
void Alquimia_PK::Setup(const Teuchos::Ptr<State>& S)
{
  Chemistry_PK::Setup(S);
 
  // Set up auxiliary chemistry data using the ChemistryEngine.
  chem_engine_->GetAuxiliaryOutputNames(aux_names_);

  for (size_t i = 0; i < aux_names_.size(); ++i) {
    std::vector<std::vector<std::string> > subname(1);
    subname[0].push_back("0");
    aux_names_[i] = Keys::getKey(domain_, aux_names_[i]);

    if (!S->HasField(aux_names_[i])) {
      S->RequireField(aux_names_[i], passwd_, subname)
       ->SetMesh(mesh_)->SetGhosted(false)
       ->SetComponent("cell", AmanziMesh::CELL, 1);
    }
  }

  if (cp_list_->isParameter("auxiliary data")) {
    auto names = cp_list_->get<Teuchos::Array<std::string> >("auxiliary data");  
    
    for (auto it = names.begin(); it != names.end(); ++it) {
      Key aux_field_name = Keys::getKey(domain_, *it);
      if (!S->HasField(aux_field_name)) {
        std::vector<std::vector<std::string> > subname(1);
        subname[0].push_back("0");
        S->RequireField(aux_field_name, passwd_, subname)
          ->SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, 1);
      }
    }
  }

  // Setup more auxiliary data
  if (!S->HasField(alquimia_aux_data_key_)) {
    int num_aux_data = chem_engine_->Sizes().num_aux_integers + chem_engine_->Sizes().num_aux_doubles;
    S->RequireField(alquimia_aux_data_key_, passwd_)
      ->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, num_aux_data);

    S->GetField(alquimia_aux_data_key_, passwd_)->set_io_vis(false);
  }
}


/* *******************************************************************
* Initialization
******************************************************************* */
void Alquimia_PK::Initialize(const Teuchos::Ptr<State>& S) 
{ 
  // initilaization using the base class
  Chemistry_PK::Initialize(S);

  if (!aux_names_.empty()) {
    aux_output_ = Teuchos::rcp(new Epetra_MultiVector(mesh_->cell_map(false), aux_names_.size()));
  } else {
    aux_output_ = Teuchos::null;
  }

  // Read XML parameters from our input file.
  XMLParameters();

  // initialize fields as soon as possible
  for (size_t i = 0; i < aux_names_.size(); ++i) {
    InitializeField(S, passwd_, aux_names_[i], 0.0);
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
  if (fabs(initial_conditions_time_ - S->time()) < 1e-8 * (1.0 + fabs(S->time()))) {
    for (auto it = chem_initial_conditions_.begin(); it != chem_initial_conditions_.end(); ++it) {
      std::string region = it->first;
      std::string condition = it->second;

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "enforcing geochemical condition \"" << condition 
                   << "\" in region \"" << region << "\"\n";
      }

      // Get the cells that belong to this region.
      int num_cells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
      AmanziMesh::Entity_ID_List cell_indices;
      mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &cell_indices);
  
      // Ensure dependencies are filled
      S_->GetFieldEvaluator(poro_key_)->HasFieldChanged(S_.ptr(), name_);
      S_->GetFieldEvaluator(fluid_den_key_)->HasFieldChanged(S_.ptr(), name_);
      S_->GetFieldEvaluator(saturation_key_)->HasFieldChanged(S_.ptr(), name_);

      // Loop over the cells.
      for (int i = 0; i < num_cells; ++i) {
        int cell = cell_indices[i];
        ierr = InitializeSingleCell(cell, condition);
      }
    }
  }

  // figure out if any of the processes threw an error, if so all processes will re-throw
  int recv = 0;
  mesh_->get_comm()->MaxAll(&ierr, &recv, 1);
  if (recv != 0) {
    Errors::Message msg("Error in Alquimia_PK::Initialize()");
    Exceptions::amanzi_throw(msg); 
  }

  // now publish auxiliary data to state
  if (aux_output_ != Teuchos::null) {
    int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

    for (int i = 0; i < aux_output_->NumVectors(); ++i) {
      auto& aux_state = *S->GetFieldData(aux_names_[i], passwd_)->ViewComponent("cell");
      for (int c = 0; c < ncells_owned; ++c) {
        aux_state[0][c] = (*aux_output_)[i][c];
      }
    }
  }

  chem_initialized_ = true;
  num_iterations_ = 0;
  num_successful_steps_ = 0;

  // verbose message
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("green") << "Initialization of PK was successful, T="
        << S->time() << vo_->reset() << std::endl << std::endl;
  }

  // S->WriteStatistics(vo_);
}


/* *******************************************************************
* This helper performs initialization on a single cell within Amanzi's state.
* It returns an error code that indicates success (0) or failure (1).
******************************************************************* */
int Alquimia_PK::InitializeSingleCell(int cell, const std::string& condition) 
{
  Teuchos::RCP<Epetra_MultiVector> tcc =
      S_->GetFieldData(tcc_key_, passwd_)->ViewComponent("cell", true);

  CopyToAlquimia(cell, tcc, alq_mat_props_, alq_state_, alq_aux_data_);

  chem_engine_->EnforceCondition(condition, current_time_, alq_mat_props_, 
                                 alq_state_, alq_aux_data_, alq_aux_output_);

  CopyAlquimiaStateToAmanzi(cell, alq_mat_props_, alq_state_, alq_aux_data_,
                            alq_aux_output_, tcc);

  return 0;
}


/* *******************************************************************
* Initialization helper functions
******************************************************************* */
void Alquimia_PK::ParseChemicalConditionRegions(const Teuchos::ParameterList& param_list,
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
    if (!cond_sublist.isType<Teuchos::Array<std::string> >("regions")) {
      msg << "Alquimia_PK::ParseChemicalConditionRegions(): \n";
      msg << "  Geochemical condition '" << cond_name << "' has no valid 'regions' entry.\n";
      Exceptions::amanzi_throw(msg);
    }

    Teuchos::Array<std::string> regions = cond_sublist.get<Teuchos::Array<std::string> >("regions");
    for (size_t r = 0; r < regions.size(); ++r) {
      // We allow for cell-based and face-based regions to accommodate both 
      // initial and boundary conditions.
      if (!mesh_->valid_set_name(regions[r], AmanziMesh::CELL) &&
          !mesh_->valid_set_name(regions[r], AmanziMesh::FACE)) {
        msg << "Alquimia_PK::ParseChemicalConditionRegions(): \n";
        msg << "  Invalid region '" << regions[r] << "' given for geochemical condition '" << cond_name << "'.\n";
        Exceptions::amanzi_throw(msg);
      }
      conditions[regions[r]] = cond_name;
    }
  }
}
 

/* *******************************************************************
*
******************************************************************* */
void Alquimia_PK::XMLParameters() 
{
  Errors::Message msg;
  Teuchos::OSTab tab = vo_->getOSTab();

  // Add any geochemical conditions we find in the Chemistry section of the file.
  if (cp_list_->isParameter("geochemical conditions")) {
    Teuchos::ParameterList conditions = cp_list_->sublist("geochemical conditions");
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
        static const char* valid_types[] = {"total_aqueous", "total_sorb", "free",
                                            "mineral", "gas", "pH", "charge"};
        static int num_valid_types = 7;
        std::string type;
        for (int i = 0; i < num_valid_types; ++i) {
          if (aqueous_constraint.isParameter(valid_types[i]))
              type = std::string(valid_types[i]);
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
  if (cp_list_->isSublist("initial conditions")) {
    // ATS-style input spec -- initial conditions in the PK
    initial_conditions = Teuchos::sublist(cp_list_, "initial conditions");
  } else {
    // Amanzi-style input spec -- initial conditions in State
    initial_conditions = Teuchos::sublist(Teuchos::sublist(glist_, "state"), "initial conditions");
  }
  if (initial_conditions->isSublist("geochemical conditions")) {
    Teuchos::ParameterList& geochem_conditions = initial_conditions->sublist("geochemical conditions");
    ParseChemicalConditionRegions(geochem_conditions, chem_initial_conditions_);
    if (chem_initial_conditions_.empty()) {
      if (cp_list_->isSublist("initial conditions")) {
        msg << "Alquimia_PK::XMLParameters(): No geochemical conditions were found in \"PK->initial conditions->geochemical conditions\"";
      } else {
        msg << "Alquimia_PK::XMLParameters(): No geochemical conditions were found in \"State->initial conditions->geochemical conditions\"";
      }
      Exceptions::amanzi_throw(msg);
    }
  }

  // Other settings.
  max_time_step_ = cp_list_->get<double>("max time step (s)", 9.9e9);
  min_time_step_ = cp_list_->get<double>("min time step (s)", 9.9e9);
  prev_time_step_ = cp_list_->get<double>("initial time step (s)", std::min(min_time_step_, max_time_step_));

  if (prev_time_step_ > max_time_step_) {
    msg << "Alquimia_PK::XMLParameters(): \n";
    msg << "  Initial Time Step exceeds Max Time Step!\n";
    Exceptions::amanzi_throw(msg);
  }

  time_step_ = prev_time_step_;
  time_step_control_method_ = cp_list_->get<std::string>("time step control method", "fixed");
  num_iterations_for_time_step_cut_ = cp_list_->get<int>("time step cut threshold", 8);
  if (num_iterations_for_time_step_cut_ <= 0) {
    msg << "Alquimia_PK::XMLParameters(): \n";
    msg << "  Invalid \"time step cut threshold\": " << num_iterations_for_time_step_cut_ << " (must be > 1).\n";
    Exceptions::amanzi_throw(msg);
  }
  num_steps_before_time_step_increase_ = cp_list_->get<int>("time step increase threshold", 4);
  if (num_steps_before_time_step_increase_ <= 0) {
    msg << "Alquimia_PK::XMLParameters(): \n";
    msg << "  Invalid \"time step increase threshold\": " << num_steps_before_time_step_increase_ << " (must be > 1).\n";
    Exceptions::amanzi_throw(msg);
  }
  time_step_cut_factor_ = cp_list_->get<double>("time step cut factor", 2.0);
  if (time_step_cut_factor_ <= 1.0) {
    msg << "Alquimia_PK::XMLParameters(): \n";
    msg << "  Invalid \"time step cut factor\": " << time_step_cut_factor_ << " (must be > 1).\n";
    Exceptions::amanzi_throw(msg);
  }
  time_step_increase_factor_ = cp_list_->get<double>("time step increase factor", 1.2);
  if (time_step_increase_factor_ <= 1.0) {
    msg << "Alquimia_PK::XMLParameters(): \n";
    msg << "  Invalid \"time step increase factor\": " << time_step_increase_factor_ << " (must be > 1).\n";
    Exceptions::amanzi_throw(msg);
  }
}


/* *******************************************************************
*
******************************************************************* */
void Alquimia_PK::CopyToAlquimia(int cell,
                                 AlquimiaProperties& mat_props,
                                 AlquimiaState& state,
                                 AlquimiaAuxiliaryData& aux_data)
{
  Teuchos::RCP<const Epetra_MultiVector> tcc =
      S_->GetFieldData(tcc_key_)->ViewComponent("cell", true);
  CopyToAlquimia(cell, tcc, mat_props, state, aux_data);
}


/* *******************************************************************
*
******************************************************************* */
void Alquimia_PK::CopyToAlquimia(int cell,
                                 Teuchos::RCP<const Epetra_MultiVector> aqueous_components,
                                 AlquimiaProperties& mat_props,
                                 AlquimiaState& state,
                                 AlquimiaAuxiliaryData& aux_data)
{
  const Epetra_MultiVector& porosity = *S_->GetFieldData(poro_key_)->ViewComponent("cell", true);
  const Epetra_MultiVector& fluid_density = *S_->GetFieldData(fluid_den_key_)->ViewComponent("cell", true);

  state.water_density = fluid_density[0][cell]; 
  state.porosity = porosity[0][cell];

  for (int i = 0; i < number_aqueous_components_; i++) {
    state.total_mobile.data[i] = (*aqueous_components)[i][cell];
    if (using_sorption_) {
      const Epetra_MultiVector& sorbed = *S_->GetFieldData(total_sorbed_key_)->ViewComponent("cell");
      state.total_immobile.data[i] = sorbed[i][cell];
    } 
  }

  // minerals
  assert(state.mineral_volume_fraction.size == number_minerals_);
  assert(state.mineral_specific_surface_area.size == number_minerals_);
  assert(mat_props.mineral_rate_cnst.size == number_minerals_);

  if (number_minerals_ > 0) {
    const Epetra_MultiVector& mineral_vf = *S_->GetFieldData(min_vol_frac_key_)->ViewComponent("cell");
    const Epetra_MultiVector& mineral_ssa = *S_->GetFieldData(min_ssa_key_)->ViewComponent("cell");
    const Epetra_MultiVector& mineral_rate = *S_->GetFieldData(mineral_rate_constant_key_)->ViewComponent("cell");
    for (unsigned int i = 0; i < number_minerals_; ++i) {
      state.mineral_volume_fraction.data[i] = mineral_vf[i][cell];
      mat_props.mineral_rate_cnst.data[i] = mineral_rate[i][cell];
      state.mineral_specific_surface_area.data[i] = mineral_ssa[i][cell];
    }
  }

  // ion exchange
  assert(state.cation_exchange_capacity.size == number_ion_exchange_sites_);
  if (number_ion_exchange_sites_ > 0) {
    const Epetra_MultiVector& ion_exchange = *S_->GetFieldData(ion_exchange_sites_key_)->ViewComponent("cell");
    for (int i = 0; i < number_ion_exchange_sites_; i++) {
      state.cation_exchange_capacity.data[i] = ion_exchange[i][cell];
    }
  }
  
  // surface complexation
  if (number_sorption_sites_ > 0) {
    const Epetra_MultiVector& sorption_sites = *S_->GetFieldData(sorp_sites_key_)->ViewComponent("cell");

    assert(number_sorption_sites_ == state.surface_site_density.size);
    for (int i = 0; i < number_sorption_sites_; ++i) {
      // FIXME: Need site density names, too?
      state.surface_site_density.data[i] = sorption_sites[i][cell];
      // TODO(bandre): need to save surface complexation free site conc here!
    }
  }

  // Auxiliary data -- block copy.
  if (S_->HasField(alquimia_aux_data_key_)) {
    aux_data_ = S_->GetFieldData(alquimia_aux_data_key_, passwd_)->ViewComponent("cell");
    int num_aux_ints = chem_engine_->Sizes().num_aux_integers;
    int num_aux_doubles = chem_engine_->Sizes().num_aux_doubles;

    for (int i = 0; i < num_aux_ints; i++) {
      double* cell_aux_ints = (*aux_data_)[i];
      aux_data.aux_ints.data[i] = (int)cell_aux_ints[cell];
    }
    for (int i = 0; i < num_aux_doubles; i++) {
      double* cell_aux_doubles = (*aux_data_)[i + num_aux_ints];
      aux_data.aux_doubles.data[i] = cell_aux_doubles[cell];
    }
  }

  const Epetra_MultiVector& water_saturation = *S_->GetFieldData(saturation_key_)->ViewComponent("cell", true);

  mat_props.volume = mesh_->cell_volume(cell);
  mat_props.saturation = water_saturation[0][cell];

  // sorption isotherms
  if (using_sorption_isotherms_) {
    const Epetra_MultiVector& isotherm_kd = *S_->GetFieldData(isotherm_kd_key_)->ViewComponent("cell");
    const Epetra_MultiVector& isotherm_freundlich_n = *S_->GetFieldData(isotherm_freundlich_n_key_)->ViewComponent("cell");
    const Epetra_MultiVector& isotherm_langmuir_b = *S_->GetFieldData(isotherm_langmuir_b_key_)->ViewComponent("cell");

    for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
      mat_props.isotherm_kd.data[i] = isotherm_kd[i][cell];
      mat_props.freundlich_n.data[i] = isotherm_freundlich_n[i][cell];
      mat_props.langmuir_b.data[i] = isotherm_langmuir_b[i][cell];
    }
  }
  
  // first order reaction rate cnst
  if (number_aqueous_kinetics_ > 0) {
    const Epetra_MultiVector& aqueous_kinetics_rate = *S_->GetFieldData(first_order_decay_constant_key_)->ViewComponent("cell");
    for (unsigned int i = 0; i < number_aqueous_kinetics_; ++i) {
      mat_props.aqueous_kinetic_rate_cnst.data[i] = aqueous_kinetics_rate[i][cell];
    }
  }
}


/* *******************************************************************
*
******************************************************************* */
void Alquimia_PK::CopyAlquimiaStateToAmanzi(
    const int cell,
    const AlquimiaProperties& mat_props,
    const AlquimiaState& state,
    const AlquimiaAuxiliaryData& aux_data,
    const AlquimiaAuxiliaryOutputData& aux_output,
    Teuchos::RCP<Epetra_MultiVector> aqueous_components)
{
  CopyFromAlquimia(cell, mat_props, state, aux_data, aux_output, 
                   aqueous_components);

  // Auxiliary output.
  std::string full_name;

  if (aux_output_ != Teuchos::null) {
    int numAqueousComplexes = chem_engine_->NumAqueousComplexes();

    for (int n = 0; n < map_[0].size(); ++n) 
      (*aux_output_)[map_[0][n]][cell] = aux_output.pH;

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
x * 
******************************************************************* */
void Alquimia_PK::CopyFromAlquimia(const int cell,
                                   const AlquimiaProperties& mat_props,
                                   const AlquimiaState& state,
                                   const AlquimiaAuxiliaryData& aux_data,
                                   const AlquimiaAuxiliaryOutputData& aux_output,
                                   Teuchos::RCP<Epetra_MultiVector> aqueous_components)
{
  // If the chemistry has modified the porosity and/or density, it needs to 
  // be updated here.
  // (this->water_density())[cell] = state.water_density;
  // (this->porosity())[cell] = state.porosity;

  for (int i = 0; i < number_aqueous_components_; ++i) {
    (*aqueous_components)[i][cell] = state.total_mobile.data[i];

    if (using_sorption_) {
      const Epetra_MultiVector& sorbed = *S_->GetFieldData(total_sorbed_key_)->ViewComponent("cell");
      sorbed[i][cell] = state.total_immobile.data[i];
    }
  }

  // Free ion species.
  const Epetra_MultiVector& free_ion = *S_->GetFieldData(free_ion_species_key_)->ViewComponent("cell");
  for (int i = 0; i < number_aqueous_components_; ++i) {
    free_ion[i][cell] = aux_output.primary_free_ion_concentration.data[i];
  }

  // Mineral properties.
  if (number_minerals_ > 0) {
    const Epetra_MultiVector& mineral_vf = *S_->GetFieldData(min_vol_frac_key_)->ViewComponent("cell");
    const Epetra_MultiVector& mineral_ssa = *S_->GetFieldData(min_ssa_key_)->ViewComponent("cell");
    const Epetra_MultiVector& mineral_rate = *S_->GetFieldData( mineral_rate_constant_key_)->ViewComponent("cell");

    for (int i = 0; i < number_minerals_; ++i) {
      mineral_vf[i][cell] = state.mineral_volume_fraction.data[i];
      mineral_ssa[i][cell] = state.mineral_specific_surface_area.data[i];
      mineral_rate[i][cell] = mat_props.mineral_rate_cnst.data[i];
    }
  }

  // ion exchange
  if (number_ion_exchange_sites_ > 0) {
    const Epetra_MultiVector& ion_exchange = *S_->GetFieldData(ion_exchange_sites_key_)->ViewComponent("cell");
    for (unsigned int i = 0; i < number_ion_exchange_sites_; i++) {
      ion_exchange[i][cell] = state.cation_exchange_capacity.data[i];
    }
  }

  // surface complexation
  if (number_sorption_sites_ > 0) {
    const Epetra_MultiVector& sorption_sites = *S_->GetFieldData(sorp_sites_key_)->ViewComponent("cell");

    for (unsigned int i = 0; i < number_sorption_sites_; i++) {
      sorption_sites[i][cell] = state.surface_site_density.data[i];
    }
  }

  if (S_->HasField(alquimia_aux_data_key_)) {
    aux_data_ = S_->GetFieldData(alquimia_aux_data_key_, passwd_)->ViewComponent("cell");

    int num_aux_ints = chem_engine_->Sizes().num_aux_integers;
    int num_aux_doubles = chem_engine_->Sizes().num_aux_doubles;

    for (int i = 0; i < num_aux_ints; i++) {
      double* cell_aux_ints = (*aux_data_)[i];
      cell_aux_ints[cell] = (double)aux_data.aux_ints.data[i];
    }
    for (int i = 0; i < num_aux_doubles; i++) {
      double* cell_aux_doubles = (*aux_data_)[i + num_aux_ints];
      cell_aux_doubles[cell] = aux_data.aux_doubles.data[i];
    } 
  }

  if (using_sorption_isotherms_) {
    const Epetra_MultiVector& isotherm_kd = *S_->GetFieldData(isotherm_kd_key_)->ViewComponent("cell");
    const Epetra_MultiVector& isotherm_freundlich_n = *S_->GetFieldData(isotherm_freundlich_n_key_)->ViewComponent("cell");
    const Epetra_MultiVector& isotherm_langmuir_b = *S_->GetFieldData(isotherm_langmuir_b_key_)->ViewComponent("cell");

    for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
      isotherm_kd[i][cell] = mat_props.isotherm_kd.data[i];
      isotherm_freundlich_n[i][cell] = mat_props.freundlich_n.data[i];
      isotherm_langmuir_b[i][cell] = mat_props.langmuir_b.data[i];
    }
  }
}


/* *******************************************************************
* This helper advances the solution on a single cell within Amanzi's state.
* It returns the number of iterations taken to obtain the advanced solution, 
* or -1 if an error occurred.
******************************************************************* */
int Alquimia_PK::AdvanceSingleCell(
    double dt, Teuchos::RCP<Epetra_MultiVector>& aqueous_components,
    int cell)
{
  // Copy the state and property information from Amanzi's state within 
  // this cell to Alquimia.
  CopyToAlquimia(cell, aqueous_components, 
                 alq_mat_props_, alq_state_, alq_aux_data_);

  // Do the reaction.
  int num_iterations;
  bool success = chem_engine_->Advance(dt, alq_mat_props_, alq_state_, 
                                       alq_aux_data_, alq_aux_output_, num_iterations);
  if (not success) 
    return -1;

  // Move the information back into Amanzi's state, updating the given total concentration vector.
  CopyAlquimiaStateToAmanzi(cell, 
                            alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_,
                            aqueous_components);

  return num_iterations;
}


/* *******************************************************************
* This function advances concentrations in the auxialiry vector 
* aqueous_components_ (defined in the base class). This vector must be
* set up using routine set_aqueous_components(). Tipically, it
* contains values advected by the transport PK.
******************************************************************* */
bool Alquimia_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool failed(false);

  double dt = t_new - t_old;
  current_time_ = saved_time_ + dt;

  // If we are given a dt that is less than the one we wanted, we don't record it.
  if (dt < time_step_) {
    prev_time_step_ = time_step_;
  } else {
    prev_time_step_ = dt;
  }

  // Get the number of owned (non-ghost) cells for the mesh.
  int num_cells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  
  int max_itrs(0), avg_itrs(0), min_itrs(1000), imax(-1);

  // Ensure dependencies are filled
  S_->GetFieldEvaluator(poro_key_)->HasFieldChanged(S_.ptr(), name_);
  S_->GetFieldEvaluator(fluid_den_key_)->HasFieldChanged(S_.ptr(), name_);
  S_->GetFieldEvaluator(saturation_key_)->HasFieldChanged(S_.ptr(), name_);

  // Now loop through all the cells and advance the chemistry.
  int convergence_failure = 0;
  for (int cell = 0; cell < num_cells; ++cell) {
    int num_itrs = AdvanceSingleCell(dt, aqueous_components_, cell);
    if (num_itrs >= 0) {
      if (max_itrs < num_itrs) {
        max_itrs = num_itrs;
        imax = cell;
      }
      min_itrs = std::min(min_itrs, num_itrs);
      avg_itrs += num_itrs;
    } else {
      // Convergence failure. Compute the next time step size.
      convergence_failure = 1;
      break;
    }
  }

  // Check for convergence failure and broadcast if needed. Also agree on the maximum number 
  // of Newton iterations and its location.
  int send[3], recv[3];
  send[0] = convergence_failure;
  send[1] = max_itrs;
  send[2] = mesh_->cell_map(false).GID(imax);
  mesh_->get_comm()->MaxAll(send, recv, 3);

  int tmp(min_itrs);
  mesh_->get_comm()->MinAll(&tmp, &min_itrs, 1);

  tmp = avg_itrs;
  mesh_->get_comm()->SumAll(&tmp, &avg_itrs, 1);
  avg_itrs /= mesh_->cell_map(false).NumGlobalElements();

  if (recv[0] != 0) 
    num_successful_steps_ = 0;
  else
    num_successful_steps_++;
  num_iterations_ = recv[1];
  imax = recv[2];

  // Compute the next time step.
  ComputeNextTimeStep();

  if (recv[0] != 0) {
    Errors::Message msg;
    msg << "Failure in Alquimia_PK::AdvanceStep";
    Exceptions::amanzi_throw(msg); 
  }
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "min/avg/max Newton: " << min_itrs << "/" << avg_itrs << "/" << num_iterations_
               << ", the maximum is in cell " << imax << std::endl;
  }

  // now publish auxiliary data to state
  if (aux_output_ != Teuchos::null) {
    int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

    for (int i = 0; i < aux_output_->NumVectors(); ++i) {
      Key full_name = aux_names_[i];
      Epetra_MultiVector& aux_state = *S_->GetFieldData(full_name, passwd_)->ViewComponent("cell");
      for (int c = 0; c < ncells_owned; ++c) {
        aux_state[0][c] = (*aux_output_)[i][c];
      }
    }
  }

  return failed;
}


/* *******************************************************************
* Time step calculation based on control parameters.
******************************************************************* */
void Alquimia_PK::ComputeNextTimeStep()
{
  if (time_step_control_method_ == "simple") {
    if ((num_successful_steps_ == 0) || (num_iterations_ >= num_iterations_for_time_step_cut_)) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "Number of Newton iterations exceeds threshold (" << num_iterations_for_time_step_cut_ 
                   << ") for time step cut, cutting dT by " << time_step_cut_factor_ << std::endl;
      }
      time_step_ = prev_time_step_ / time_step_cut_factor_;
    }
    else if (num_successful_steps_ >= num_steps_before_time_step_increase_) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "Number of successful steps exceeds threshold (" << num_steps_before_time_step_increase_ 
                   << ") for time step increase, growing dT by " << time_step_increase_factor_ << std::endl;
      }
      time_step_ = prev_time_step_ * time_step_increase_factor_;
      num_successful_steps_ = 0;
    }
  }

  if (time_step_ > max_time_step_)
    time_step_ = max_time_step_;
  /* else if (time_step_ > min_time_step_)
     time_step_ = min_time_step_; */
}


/* *******************************************************************
* TBW
******************************************************************* */
void Alquimia_PK::CopyFieldstoNewState(const Teuchos::RCP<State>& S_next)
{
  Chemistry_PK::CopyFieldstoNewState(S_next);

  std::vector<std::string> aux_names;
  chem_engine_->GetAuxiliaryOutputNames(aux_names);

  for (size_t i = 0; i < aux_names.size(); ++i) {
    std::vector<std::vector<std::string> > subname(1);
    subname[0].push_back("0");
    aux_names[i] = Keys::getKey(domain_, aux_names[i]);
    if (S_->HasField(aux_names[i])&&S_next->HasField(aux_names[i])) {
      *S_next->GetFieldData(aux_names[i], passwd_)->ViewComponent("cell", false) =
        *S_->GetFieldData(aux_names[i], passwd_)->ViewComponent("cell", false);
    }
  }

  if (cp_list_->isParameter("auxiliary data")) {
    Teuchos::Array<std::string> names = 
      cp_list_->get<Teuchos::Array<std::string> >("auxiliary data");  
    
    for (auto it = names.begin(); it != names.end(); ++it) {
      Key aux_field_name = Keys::getKey(domain_, *it);
      if (S_->HasField(aux_field_name) && S_next->HasField(aux_field_name)) {
        *S_next->GetFieldData(aux_field_name, passwd_)->ViewComponent("cell", false) =
          *S_->GetFieldData(aux_field_name, passwd_)->ViewComponent("cell", false);
      }
    }
  }

  if (S_->HasField(alquimia_aux_data_key_) && S_next->HasField(alquimia_aux_data_key_)) {
    *S_next->GetFieldData(alquimia_aux_data_key_, passwd_)->ViewComponent("cell", false) =
      *S_->GetFieldData(alquimia_aux_data_key_, passwd_)->ViewComponent("cell", false);
  }
}


/* *******************************************************************
* The MPC will call this function to signal to the process kernel that 
* it has accepted the state update, thus, the PK should update
* possible auxilary state variables here
******************************************************************* */
void Alquimia_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) 
{
  if (S_ != S) CopyFieldstoNewState(S);
  saved_time_ = t_new;
}


/* *******************************************************************
*
******************************************************************* */
Teuchos::RCP<Epetra_MultiVector> Alquimia_PK::extra_chemistry_output_data() 
{
  // This vector is updated during the initialization and advance of 
  // the geochemistry, so we simply return it here.
  return aux_output_;
}


/* *******************************************************************
* Auxiliary map from beacon to aux_names
******************************************************************* */
void Alquimia_PK::InitializeAuxNamesMap_()
{
  map_.resize(7);

  std::string full_name;
  int numAqueousComplexes = chem_engine_->NumAqueousComplexes();

  for (int i = 0; i < aux_names_.size(); i++) {
    if (aux_names_.at(i) == "pH") {
      map_[0].push_back(i);
    }
    else if (aux_names_.at(i).find("mineral_saturation_index") != std::string::npos) {
      for (int j = 0; j < mineral_names_.size(); ++j) {
        full_name = Keys::getKey(domain_,std::string("mineral_saturation_index_") + mineral_names_[j]);
        if (aux_names_.at(i) == full_name) {
          map_[1].push_back(i);
        }
      }
    }
    else if (aux_names_.at(i).find("mineral_reaction_rate") != std::string::npos) {
      for (int j = 0; j < mineral_names_.size(); ++j) {
        full_name = Keys::getKey(domain_,std::string("mineral_reaction_rate_") + mineral_names_[j]);
        if (aux_names_.at(i) == full_name) {
          map_[2].push_back(i);
        }
      }
    }
    else if (aux_names_.at(i).find("primary_free_ion_concentration") != std::string::npos) {
      for (int j = 0; j < primary_names_.size(); ++j) {
        full_name = Keys::getKey(domain_,std::string("primary_free_ion_concentration_") + primary_names_[j]);
        if (aux_names_.at(i) == full_name) {
          map_[3].push_back(i);
        }
      }
    }
    else if (aux_names_.at(i).find(primary_activity_coeff_key_) != std::string::npos) {
      for (int j = 0; j < primary_names_.size(); ++j) {
        full_name = Keys::getKey(domain_,std::string("primary_activity_coeff_") + primary_names_[j]);
        if (aux_names_.at(i) == full_name) {
          map_[4].push_back(i);
        }
      }
    }
    else if (aux_names_.at(i).find("secondary_free_ion_concentration") != std::string::npos) {
      for (int j = 0; j < numAqueousComplexes; ++j) {
        char num_str[16];
        snprintf(num_str, 15, "%d", j);
        full_name = Keys::getKey(domain_,std::string("secondary_free_ion_concentration_") + std::string(num_str));
        if (aux_names_.at(i) == full_name) {
          map_[5].push_back(i);
        }
      }
    }
    else if (aux_names_.at(i).find(secondary_activity_coeff_key_) != std::string::npos) {
      for (int j = 0; j < numAqueousComplexes; ++j) {
        char num_str[16];
        snprintf(num_str, 15, "%d", j);
        full_name = Keys::getKey(domain_,std::string("secondary_activity_coeff_") + std::string(num_str));
        if (aux_names_.at(i) == full_name) {
          map_[6].push_back(i);
        }
      }
    }
  }
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
