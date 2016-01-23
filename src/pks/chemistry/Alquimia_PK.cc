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

#include <set>
#include <string>
#include <algorithm>

// TPLs
#include "boost/mpi.hpp"
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
Alquimia_PK::Alquimia_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                         Teuchos::RCP<Chemistry_State> chem_state,
                         Teuchos::RCP<ChemistryEngine> chem_engine,
                         Teuchos::RCP<State> S,
                         Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
    glist_(glist),
    max_time_step_(9.9e9),
    chemistry_state_(chem_state),
    chem_initialized_(false),
    chem_engine_(chem_engine),
    current_time_(0.0),
    saved_time_(0.0) 
{
  S_ = S;
  mesh_ = mesh;

  // We need the chemistry list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  cp_list_ = Teuchos::sublist(pk_list, "Chemistry", true);

  InitializeMinerals(cp_list_);

  // verbosity object
  vo_ = Teuchos::rcp(new VerboseObject("Chemistry_PK", *cp_list_));
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
void Alquimia_PK::Setup()
{
  Chemistry_PK::Setup();
}


/* *******************************************************************
* This helper performs initialization on a single cell within Amanzi's state.
* It returns an error code that indicates success (0) or failure (1).
******************************************************************* */
int Alquimia_PK::InitializeSingleCell(int cell_index, const std::string& condition) 
{
  Teuchos::RCP<Epetra_MultiVector> tcc =
      S_->GetFieldData("total_component_concentration", passwd_)->ViewComponent("cell");

  // Copy the state and material information from Amanzi's state within 
  // this cell to Alquimia.
  CopyAmanziStateToAlquimia(cell_index, tcc,
                            alq_mat_props_, alq_state_, alq_aux_data_);

  // Do the initialization.
  chem_engine_->EnforceCondition(condition, current_time_, alq_mat_props_, 
                                 alq_state_, alq_aux_data_, alq_aux_output_);

  // Move the information back into Amanzi's state.
  CopyAlquimiaStateToAmanzi(cell_index,
                            alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_,
                            tcc);

  return 0;
}


/* *******************************************************************
* 
******************************************************************* */
void Alquimia_PK::InitializeChemistry() 
{
  Errors::Message msg;

  // Read XML parameters from our input file.
  XMLParameters();

  // Initialize the data structures that we will use to traffic data between 
  // Amanzi and Alquimia.
  chem_engine_->InitState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);

  // Now loop through all the regions and initialize.
  int ierr = 0;
  for (std::map<std::string, std::string>::const_iterator cond_iter = chem_initial_conditions_.begin(); 
       cond_iter != chem_initial_conditions_.end(); ++cond_iter) {
    std::string region_name = cond_iter->first;
    std::string condition = cond_iter->second;

    // Get the cells that belong to this region.
    unsigned int num_cells = mesh_->get_set_size(region_name, AmanziMesh::CELL, AmanziMesh::OWNED);
    AmanziMesh::Entity_ID_List cell_indices;
    mesh_->get_set_entities(region_name, AmanziMesh::CELL, AmanziMesh::OWNED, &cell_indices);
  
    // Loop over the cells.
    for (unsigned int i = 0; i < num_cells; ++i) {
      int cell = cell_indices[i];
      ierr = InitializeSingleCell(cell, condition);
    }
  }

  // figure out if any of the processes threw an error, if so all processes will re-throw
  int recv = 0;
  mesh_->get_comm()->MaxAll(&ierr, &recv, 1);
  if (recv != 0) {
    msg << "Error in Alquimia_PK::InitializeChemistry 1";
    Exceptions::amanzi_throw(msg); 
  }

  // now publish auxiliary data to state
  if (aux_output_ != Teuchos::null) {
    for (int i = 0; i < aux_output_->NumVectors(); ++i) {
      Epetra_MultiVector& aux_state = *S_->GetFieldData(aux_names_[i], passwd_)->ViewComponent("cell");
      aux_state[0] = (*aux_output_)[i];
    }
  }

  chem_initialized_ = true;
  num_iterations_ = 0;
  num_successful_steps_ = 0;
}


/* *******************************************************************
* Initialization helper functions
******************************************************************* */
void Alquimia_PK::ParseChemicalConditionRegions(const Teuchos::ParameterList& param_list,
                                                std::map<std::string, std::string>& conditions)
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
void Alquimia_PK::XMLParameters(void) 
{
  Errors::Message msg;

  // We retrieve the names of the auxiliary output data from the chemistry engine--we don't rely on 
  // the Auxiliary Data parameter list.
  chem_engine_->GetAuxiliaryOutputNames(aux_names_);
  if (!aux_names_.empty()) {
    aux_output_ = Teuchos::rcp(new Epetra_MultiVector(mesh_->cell_map(false), aux_names_.size()));
  }
  else {
    aux_output_ = Teuchos::null;
  }

  // Add any geochemical conditions we find in the Chemistry section of the file.
  if (cp_list_->isParameter("geochemical conditions")) {
    Teuchos::ParameterList conditions = cp_list_->sublist("geochemical conditions");
    for (Teuchos::ParameterList::ConstIterator iter = conditions.begin();
         iter != conditions.end(); ++iter) {
      // This parameter list contains sublists, each corresponding to a
      // geochemical condition. 
      std::string cond_name = conditions.name(iter);
      assert(conditions.isSublist(cond_name));
      const Teuchos::ParameterList& cond_sublist = conditions.sublist(cond_name);

      // Create the entry for this geochemical condition within the chemistry engine,
      // overwriting any previous definition.
      chem_engine_->CreateCondition(cond_name);

      // Now mine the entry for details.
      for (Teuchos::ParameterList::ConstIterator iter2 = cond_sublist.begin();
           iter2 != cond_sublist.end(); ++iter2) {
        std::string species_name = cond_sublist.name(iter2);
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
              msg << "No mineral or gas species has been given to equilibrate with species '" << species_name << "'.\n";
              Exceptions::amanzi_throw(msg);
            }
          }  
        }
        if (!type.empty()) {
          // It's a valid aqueous constraint, so we add it to the chemistry engine
          // under the current geochemical condition.
          chem_engine_->AddAqueousConstraint(cond_name, species_name, type, equilibrate_name);
        }
        else {
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
  if (!glist_->isSublist("State")) {
    msg << "Alquimia_PK::XMLParameters(): \n";
    msg << "  No 'State' sublist was found!\n";
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::ParameterList& state_list = glist_->sublist("State");
  Teuchos::ParameterList& initial_conditions = state_list.sublist("initial conditions");
  if (!initial_conditions.isSublist("geochemical conditions")) {
    msg << "Alquimia_PK::XMLParameters(): \n";
    msg << "  No 'geochemical conditions' list was found in 'State->initial conditions'!\n";
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::ParameterList geochem_conditions = initial_conditions.sublist("geochemical conditions");
  ParseChemicalConditionRegions(geochem_conditions, chem_initial_conditions_);
  if (chem_initial_conditions_.empty()) {
    msg << "Alquimia_PK::XMLParameters(): \n";
    msg << "  No geochemical conditions were found in 'State->initial conditions->geochemical conditions'!\n";
    Exceptions::amanzi_throw(msg);
  }

  // Other settings.
  max_time_step_ = cp_list_->get<double>("max time step (s)", 9.9e9);
  min_time_step_ = cp_list_->get<double>("min time step (s)", 9.9e9);
  prev_time_step_ = cp_list_->get<double>("initial time step (s)", std::min(min_time_step_, max_time_step_));
  /*
  if ((min_time_step_ == 9.9e9) && (prev_time_step_ < 9.9e9))
    min_time_step_ = prev_time_step_;
  else if ((min_time_step_ == 9.9e9) && (max_time_step_ < 9.9e9))
    min_time_step_ = max_time_step_;
  else if (min_time_step_ > max_time_step_) {
    msg << "Alquimia_PK::XMLParameters(): \n";
    msg << "  Min Time Step exceeds Max Time Step!\n";
    Exceptions::amanzi_throw(msg);
  }
  if ((min_time_step_ < max_time_step_) && (prev_time_step_ == max_time_step_))
    prev_time_step_ = min_time_step_;
  */
  if (prev_time_step_ > max_time_step_) {
    msg << "Alquimia_PK::XMLParameters(): \n";
    msg << "  Initial Time Step exceeds Max Time Step!\n";
    Exceptions::amanzi_throw(msg);
  }
  /*
   if (prev_time_step_ < min_time_step_) {
    msg << "Alquimia_PK::XMLParameters(): \n";
    msg << "  Initial Time Step is smaller than Min Time Step!\n";
    Exceptions::amanzi_throw(msg);
    }
  */

  time_step_ = prev_time_step_;
  time_step_control_method_ = cp_list_->get<std::string>("time step control method", "fixed");
  num_iterations_for_time_step_cut_ = cp_list_->get<int>("time step cut threshold", 8);
  if (num_iterations_for_time_step_cut_ <= 0)
  {
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
void Alquimia_PK::CopyAmanziStateToAlquimia(const int cell_id,
                                            Teuchos::RCP<const Epetra_MultiVector> aqueous_components,
                                            AlquimiaMaterialProperties& mat_props,
                                            AlquimiaState& state,
                                            AlquimiaAuxiliaryData& aux_data)
{
  chemistry_state_->CopyToAlquimia(cell_id, aqueous_components, mat_props, state, aux_data);
}


/* *******************************************************************
*
******************************************************************* */
void Alquimia_PK::CopyAlquimiaStateToAmanzi(
    const int cell_id,
    const AlquimiaMaterialProperties& mat_props,
    const AlquimiaState& state,
    const AlquimiaAuxiliaryData& aux_data,
    const AlquimiaAuxiliaryOutputData& aux_output,
    Teuchos::RCP<Epetra_MultiVector> total_component_concentration)
{
  chemistry_state_->CopyFromAlquimia(cell_id, mat_props, state, aux_data, aux_output, 
                                     total_component_concentration);

  // Auxiliary output.
  if (aux_output_ != Teuchos::null) {
    std::vector<std::string> mineralNames, primaryNames;
    chem_engine_->GetMineralNames(mineralNames);
    chem_engine_->GetPrimarySpeciesNames(primaryNames);
    int numAqueousComplexes = chem_engine_->NumAqueousComplexes();
    for (unsigned int i = 0; i < aux_names_.size(); i++) {
      if (aux_names_.at(i) == "pH") {
        double* cell_aux_output = (*aux_output_)[i];
        cell_aux_output[cell_id] = aux_output.pH;
      }
      else if (aux_names_.at(i).find("mineral_saturation_index") != std::string::npos) {
        for (int j = 0; j < mineralNames.size(); ++j) {
          std::string full_name = std::string("mineral_saturation_index_") + mineralNames[j];
          if (aux_names_.at(i) == full_name) {
            double* cell_aux_output = (*aux_output_)[i];
            cell_aux_output[cell_id] = aux_output.mineral_saturation_index.data[j];
          }
        }
      }
      else if (aux_names_.at(i).find("mineral_reaction_rate") != std::string::npos) {
        for (int j = 0; j < mineralNames.size(); ++j) {
          std::string full_name = std::string("mineral_reaction_rate_") + mineralNames[j];
          if (aux_names_.at(i) == full_name) {
            double* cell_aux_output = (*aux_output_)[i];
            cell_aux_output[cell_id] = aux_output.mineral_reaction_rate.data[j];
          }
        }
      }
      else if (aux_names_.at(i) == "primary_free_ion_concentration") {
        for (int j = 0; j < primaryNames.size(); ++j) {
          std::string full_name = std::string("primary_free_ion_concentration_") + primaryNames[j];
          if (aux_names_.at(i) == full_name) {
            double* cell_aux_output = (*aux_output_)[i];
            cell_aux_output[cell_id] = aux_output.primary_free_ion_concentration.data[j];
          }
        }
      }
      else if (aux_names_.at(i) == "primary_activity_coeff") {
        for (int j = 0; j < primaryNames.size(); ++j) {
          std::string full_name = std::string("primary_activity_coeff_") + primaryNames[j];
          if (aux_names_.at(i) == full_name) {
            double* cell_aux_output = (*aux_output_)[i];
            cell_aux_output[cell_id] = aux_output.primary_activity_coeff.data[j];
          }
        }
      }
      else if (aux_names_.at(i) == "secondary_free_ion_concentration") {
        for (int j = 0; j < numAqueousComplexes; ++j) {
          char num_str[16];
          snprintf(num_str, 15, "%d", j);
          std::string full_name = std::string("secondary_free_ion_concentration_") + std::string(num_str);
          if (aux_names_.at(i) == full_name) {
            double* cell_aux_output = (*aux_output_)[i];
            cell_aux_output[cell_id] = aux_output.secondary_free_ion_concentration.data[j];
          }
        }
      }
      else if (aux_names_.at(i) == "secondary_activity_coeff") {
        for (int j = 0; j < numAqueousComplexes; ++j) {
          char num_str[16];
          snprintf(num_str, 15, "%d", j);
          std::string full_name = std::string("secondary_activity_coeff_") + std::string(num_str);
          if (aux_names_.at(i) == full_name) {
            double* cell_aux_output = (*aux_output_)[i];
            cell_aux_output[cell_id] = aux_output.secondary_activity_coeff.data[j];
          }
        }
      }
    }
  }
}


/* *******************************************************************
* This helper advances the solution on a single cell within Amanzi's state.
* It returns the number of iterations taken to obtain the advanced solution, 
* or -1 if an error occurred.
******************************************************************* */
int Alquimia_PK::AdvanceSingleCell(
    double delta_time, Teuchos::RCP<Epetra_MultiVector> total_component_concentration,
    int cell_index)
{
  // Copy the state and material information from Amanzi's state within 
  // this cell to Alquimia.
  CopyAmanziStateToAlquimia(cell_index, 
                            total_component_concentration, 
                            alq_mat_props_, alq_state_, alq_aux_data_);

  // Do the reaction.
  int num_iterations;
  bool success = chem_engine_->Advance(delta_time, alq_mat_props_, alq_state_, 
                                       alq_aux_data_, alq_aux_output_, num_iterations);
  if (not success) 
    return -1;

  // Move the information back into Amanzi's state, updating the given total concentration vector.
  CopyAlquimiaStateToAmanzi(cell_index, 
                            alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_,
                            total_component_concentration);

  return num_iterations;
}


/* *******************************************************************
* - the MPC will call this function to advance the state with this
* particular process kernel
*
* - this is how to get the total component concentration
* CS->get_total_component_concentration()
*
* - please update the argument to this function called tcc_star
* with the result of your chemistry computation which is the total
* component concentration ^star
*
* - see the Chemistry_State for the other available data in the
* chemistry specific state
******************************************************************* */
void Alquimia_PK::Advance(const double& delta_time,
                          Teuchos::RCP<Epetra_MultiVector> total_component_concentration) 
{
  Errors::Message msg;

  current_time_ = saved_time_ + delta_time;

  // If we are given a dt that is less than the one we wanted, we don't record it.
  if (delta_time < time_step_) {
    prev_time_step_ = time_step_;
  }
  else {
    prev_time_step_ = delta_time;
  }

  // Get the number of owned (non-ghost) cells for the mesh.
  unsigned int num_cells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  
  int max_iterations (0), min_iterations(10000000), ave_iterations(0);
  int imax(-1), imin(-1);

  // Now loop through all the cells and advance the chemistry.
  int ierr = 0;
  int convergence_failure = 0;
  for (int cell = 0; cell < num_cells; ++cell) {
    int num_iterations = AdvanceSingleCell(delta_time, total_component_concentration, cell);
    if (num_iterations >= 0) {
      if (max_iterations < num_iterations) {
        max_iterations = num_iterations;
        imax = cell;
      }
      if (min_iterations > num_iterations) {
        min_iterations = num_iterations;
        imin = cell;
      }
      ave_iterations += num_iterations;
    } 
    else {
      // Convergence failure. Compute the next time step size.
      convergence_failure = 1;
      break;
    }
  }

  // Check for convergence failure and broadcast if needed. Also agree on the maximum number 
  // of Newton iterations and its location.
  int send[3], recv[3];
  send[0] = convergence_failure;
  send[1] = max_iterations;
  send[2] = imax;
  mesh_->get_comm()->MaxAll(send, recv, 3);
  if (recv[0] != 0) 
    num_successful_steps_ = 0;
  else
    num_successful_steps_++;
  num_iterations_ = recv[1];
  imax = recv[2];

  // Compute the next time step.
  ComputeNextTimeStep();

  if (recv[0] != 0) {
    msg << "Failure in Alquimia_PK::Advance";
    Exceptions::amanzi_throw(msg); 
  }
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Advanced after maximum of " << num_iterations_
               << " Newton iterations in cell " << imax << "." << std::endl;
  }

  // now publish auxiliary data to state
  if (aux_output_ != Teuchos::null) {
    for (int i=0; i<aux_output_->NumVectors(); ++i) {
      Epetra_MultiVector& aux_state = *S_->GetFieldData(aux_names_[i], passwd_)->ViewComponent("cell");
      aux_state[0] = (*aux_output_)[i];
    }
  }
}


/* *******************************************************************
*
******************************************************************* */
void Alquimia_PK::ComputeNextTimeStep()
{
  if (time_step_control_method_ == "simple") {
    if ((num_successful_steps_ == 0) || (num_iterations_ >= num_iterations_for_time_step_cut_)) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "Number of Newton iterations exceeds threshold (" << num_iterations_for_time_step_cut_ << ") for time step cut, cutting dT by " << time_step_cut_factor_ << std::endl;
      }
      time_step_ = prev_time_step_ / time_step_cut_factor_;
    }
    else if (num_successful_steps_ >= num_steps_before_time_step_increase_) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "Number of successful steps exceeds threshold (" << num_steps_before_time_step_increase_ << ") for time step increase, growing dT by " << time_step_increase_factor_ << std::endl;
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
* The MPC will call this function to signal to the process kernel that 
* it has accepted the state update, thus, the PK should update
* possible auxilary state variables here
******************************************************************* */
void Alquimia_PK::CommitState(Teuchos::RCP<Chemistry_State> chem_state,
                              const double& time) 
{
  saved_time_ = time;
}


/* *******************************************************************
*
******************************************************************* */
Teuchos::RCP<Epetra_MultiVector> Alquimia_PK::get_extra_chemistry_output_data() 
{
  // This vector is updated during the initialization and advance of 
  // the geochemistry, so we simply return it here.
  return aux_output_;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
