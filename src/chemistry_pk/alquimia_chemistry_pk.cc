/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "alquimia_chemistry_pk.hh"

#include <set>
#include <string>
#include <algorithm>

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "errors.hh"
#include "exceptions.hh"

#include "boost/mpi.hpp"

namespace Amanzi {
namespace AmanziChemistry {

/*******************************************************************************
 **
 **  Purpose: Trilinos based process kernel for chemistry
 **
 **  Notes:
 **
 **    - all the actual geochemistry calculations live in the chemistry library. 
 **
 **    - The process kernel stores the instance of the chemistry object
 **    and drives the chemistry calculations on a cell by cell
 **    basis. It handles the movement of data back and forth between
 **    the amanzi memory and the chemistry library data structures.
 **
 **    - chemistry_state will always hold the state at the begining of
 **    the time step. should not (can not?) be changed by chemistry.
 **
 **    - when advance is called, total_component_concentration_star
 **    holds the value of component concentrations after transport!
 **
 **    - where do we write to when advance state is done? tcc is read
 **    only, do we want to write over the values in tcc_star?
 **
 **    - when commit_state is called, we get a new Chemistry_State
 **    object which will hold the final info for the end of the time
 **    step. We can use it if we want, to update our internal data.
 **
 **    - The State and Chemistry_State objects have
 **    total_component_concentrations in a multi vector. The data is
 **    stored such that:
 **
 **      double* foo = total_component_concetration[i]
 *
 **    where foo refers to a vector of component concentrations for a
 **    single component.
 **
 **
 ******************************************************************************/

Alquimia_Chemistry_PK::Alquimia_Chemistry_PK(const Teuchos::ParameterList& param_list,
                                             Teuchos::RCP<Chemistry_State> chem_state,
                                             Teuchos::RCP<ChemistryEngine> chem_engine)
    : max_time_step_(9.9e9),
      chemistry_state_(chem_state),
      main_param_list_(param_list),
      chem_param_list_(),
      chem_initialized_(false),
      chem_engine_(chem_engine),
      current_time_(0.0),
      saved_time_(0.0) 
{
  // We need the top-level parameter list.
  chem_param_list_ = main_param_list_.sublist("Chemistry");
}  // end Alquimia_Chemistry_PK()

Alquimia_Chemistry_PK::~Alquimia_Chemistry_PK() 
{
  // Destroy ansilary data structures.
  chem_engine_->FreeState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);
}  // end ~Alquimia_Chemistry_PK()

// This helper performs initialization on a single cell within Amanzi's state.
// It returns an error code that indicates success (0) or failure (1).
int Alquimia_Chemistry_PK::InitializeSingleCell(int cellIndex, const std::string& condition) 
{
  // Copy the state and material information from Amanzi's state within 
  // this cell to Alquimia.
  CopyAmanziStateToAlquimia(cellIndex, 
                            chemistry_state_->total_component_concentration(), 
                            alq_mat_props_, alq_state_, alq_aux_data_);

  // Do the initialization.
  chem_engine_->EnforceCondition(condition, current_time_, alq_mat_props_, 
                                 alq_state_, alq_aux_data_, alq_aux_output_);

  // Move the information back into Amanzi's state.
  CopyAlquimiaStateToAmanzi(cellIndex, alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_,
                            chemistry_state_->total_component_concentration());

  return 0;
}

void Alquimia_Chemistry_PK::InitializeChemistry(void) 
{
  Errors::Message msg;

  // Read XML parameters from our input file.
  XMLParameters();

  // Initialize the data structures that we will use to traffic data between 
  // Amanzi and Alquimia.
  chem_engine_->InitState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);

  // Now loop through all the regions and initialize.
  int ierr = 0;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = chemistry_state_->mesh_maps();
  for (std::map<std::string, std::string>::const_iterator cond_iter = chem_initial_conditions_.begin(); 
       cond_iter != chem_initial_conditions_.end(); ++cond_iter)
  {
    std::string region_name = cond_iter->first;
    std::string condition = cond_iter->second;

    // Get the cells that belong to this region.
    unsigned int num_cells = mesh->get_set_size(region_name, AmanziMesh::CELL, AmanziMesh::OWNED);
    AmanziMesh::Entity_ID_List cell_indices;
    mesh->get_set_entities(region_name, AmanziMesh::CELL, AmanziMesh::OWNED, &cell_indices);
  
    // Loop over the cells.
    for (unsigned int i = 0; i < num_cells; ++i) 
    {
      int cell = cell_indices[i];
      ierr = InitializeSingleCell(cell, condition);
    }
  }

  // figure out if any of the processes threw an error, if so all processes will re-throw
  int recv = 0;
  chemistry_state_->mesh_maps()->get_comm()->MaxAll(&ierr, &recv, 1);
  if (recv != 0) 
  {
    msg << "Error in Alquimia_Chemistry_PK::InitializeChemistry 1";
    Exceptions::amanzi_throw(msg); 
  }

  chem_initialized_ = true;
}  // end InitializeChemistry()

/*******************************************************************************
 **
 **  initialization helper functions
 **
 ******************************************************************************/
void Alquimia_Chemistry_PK::ParseChemicalConditionRegions(const Teuchos::ParameterList& param_list,
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
    if (!cond_sublist.isType<Teuchos::Array<std::string> >("regions"))
    {
      msg << "Alquimia_Chemistry_PK::ParseChemicalConditionRegions(): \n";
      msg << "  Geochemical condition '" << cond_name << "' has no valid 'regions' entry.\n";
      Exceptions::amanzi_throw(msg);
    }
    Teuchos::Array<std::string> regions = cond_sublist.get<Teuchos::Array<std::string> >("regions");
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = chemistry_state_->mesh_maps();
    for (size_t r = 0; r < regions.size(); ++r)
    {
      // We allow for cell-based and face-based regions to accommodate both 
      // initial and boundary conditions.
      if (!mesh->valid_set_name(regions[r], AmanziMesh::CELL) &&
          !mesh->valid_set_name(regions[r], AmanziMesh::FACE))
      {
        msg << "Alquimia_Chemistry_PK::ParseChemicalConditionRegions(): \n";
        msg << "  Invalid region '" << regions[r] << "' given for geochemical condition '" << cond_name << "'.\n";
        Exceptions::amanzi_throw(msg);
      }
      conditions[regions[r]] = cond_name;
    }
  }
}

void Alquimia_Chemistry_PK::XMLParameters(void) 
{
  Errors::Message msg;

  // NOTE that our parameter list should be the top-level parameter list "Main", 
  // not the "Chemistry" one used by the native Amanzi Chemistry PK.

  // Auxiliary Data
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
      else if ((*name == "mineral_saturation_index") || (*name == "mineral_reaction_rate"))
      {
        // We need one of these for each mineral.
        std::vector<std::string> mineralNames;
        chem_engine_->GetMineralNames(mineralNames);
        for (int i = 0; i < mineralNames.size(); ++i)
        {
          std::string aux_name = *name + string("_") + mineralNames[i];
          aux_names_.push_back(*name);
        }
      }
      else if ((*name == "primary_free_ion_concentration") || (*name == "primary_activity_coeff"))
      {
        // We need one of these for each primary species.
        std::vector<std::string> primaryNames;
        chem_engine_->GetPrimarySpeciesNames(primaryNames);
        for (int i = 0; i < primaryNames.size(); ++i)
        {
          std::string aux_name = *name + string("_") + primaryNames[i];
          aux_names_.push_back(*name);
        }
      }
      else if ((*name == "secondary_free_ion_concentration") || (*name == "secondary_activity_coeff"))
      {
        // We need one of these for each secondary species, which aren't named.
        for (int i = 0; i < chem_engine_->NumAqueousComplexes(); ++i)
        {
          char num_str[16];
          snprintf(num_str, 15, "%d", i);
          std::string aux_name = *name + string("_") + string(num_str);
          aux_names_.push_back(*name);
        }
        aux_names_.push_back(*name);
      } 
      else 
      {
        std::stringstream message;
        message << "Alquimia_ChemistryPK::XMLParameters(): unknown value in 'Auxiliary Data' list: " 
                << *name << std::endl;
        chem_out->Write(kWarning, message);
      }
    }
    if (names.size() > 0)
    {
      aux_output_ = Teuchos::rcp(new Epetra_MultiVector(
           chemistry_state_->mesh_maps()->cell_map(false), names.size()));
    }
    else
    {
      aux_output_ = Teuchos::null;
    }
  }

  // Add any geochemical conditions we find in the Chemistry section of the file.
  if (chem_param_list_.isParameter("Geochemical Conditions")) 
  {
    Teuchos::ParameterList conditions = chem_param_list_.sublist("Geochemical Conditions");
    for (Teuchos::ParameterList::ConstIterator iter = conditions.begin();
         iter != conditions.end(); ++iter)
    {
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
           iter2 != cond_sublist.end(); ++iter2)
      {
        std::string species_name = cond_sublist.name(iter2);
        assert(cond_sublist.isSublist(species_name));
        const Teuchos::ParameterList& aqueous_constraint = cond_sublist.sublist(species_name);
        
        // If the primary species has an associated mineral, we need to retrieve its information.
        std::string mineral_name;
        if (aqueous_constraint.isParameter("species"))
        {
          mineral_name = aqueous_constraint.get<std::string>("species");

          // The information for this mineral species must appear alongside the 
          // aqueous constraint in the geochemical condition.
          if (!cond_sublist.isSublist(mineral_name))
          {
            msg << "Aqueous constraint for " << species_name << " requires a mineral constraint for\n";
            msg << mineral_name << ", but no such constraint is present.\n";
            Exceptions::amanzi_throw(msg);
          }

          const Teuchos::ParameterList& mineral_constraint = cond_sublist.sublist(mineral_name);

          if (!mineral_constraint.isParameter("Volume Fraction") || 
              !mineral_constraint.isParameter("Specific Surface Area"))
          {
            msg << "Mineral constraint for " << mineral_name << " requires 'Volume Fraction' and\n";
            msg << "'Specific Surface Area' to be specified, but they are not.\n";
            Exceptions::amanzi_throw(msg);
          }
          double volume_fraction = mineral_constraint.get<double>("Volume Fraction");
          if (volume_fraction < 0.0)
          {
            msg << "Negative volume fraction given for mineral constraint '" << mineral_name << "'.\n";
            Exceptions::amanzi_throw(msg);
          }
          double specific_surface_area = mineral_constraint.get<double>("Specific Surface Area");
          if (specific_surface_area < 0.0)
          {
            msg << "Negative specific surface area given for mineral constraint '" << mineral_name << "'.\n";
            Exceptions::amanzi_throw(msg);
          }
          // Convert specific surface area to m**2/m**3 from cm**2/cm**3.
          specific_surface_area *= 100.0;

          // Add the mineral constraint to the chemistry engine.
          chem_engine_->AddMineralConstraint(cond_name, mineral_name, volume_fraction, specific_surface_area);
        }

        // What kind of aqueous constraint do we have on this species?
        static const char* valid_types[] = {"total_aqueous", "total_sorb", "free",
                                            "mineral", "gas", "pH", "charge"};
        static int num_valid_types = 7;
        std::string type;
        for (int i = 0; i < num_valid_types; ++i)
        {
          if (aqueous_constraint.isParameter(valid_types[i]))
            type = std::string(valid_types[i]);
        }
        if (!type.empty())
        {
          // It's a valid aqueous constraint, so we add it to the chemistry engine
          // under the current geochemical condition.
          chem_engine_->AddAqueousConstraint(cond_name, species_name, type, mineral_name);
        }
        else
        {
          // This is either a mineral constraint or an invalid aqueous constraint.
          if (aqueous_constraint.isParameter("Volume Fraction") && 
              aqueous_constraint.isParameter("Specific Surface Area"))
          {
            // It's a mineral constraint. We process these separately above.
            continue;
          }
          else
          {
            // We have an invalid aqueous contraint.
            msg << "Invalid aqueous constraint type for " << species_name << ".\n";
            msg << "Valid types are total_aqueous, total_sorb, free, mineral, gas, pH, and charge.\n";
            Exceptions::amanzi_throw(msg);
          }
        }
      }
    }
  }

  // Now associate regions with chemical conditions based on initial 
  // condition specifications in the file.
  if (!main_param_list_.isSublist("State"))
  {
    msg << "Alquimia_Chemistry_PK::XMLParameters(): \n";
    msg << "  No 'State' sublist was found!\n";
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::ParameterList state_list = main_param_list_.sublist("State");
  Teuchos::ParameterList initial_conditions = state_list.sublist("initial conditions");
  if (!initial_conditions.isSublist("geochemical conditions"))
  {
    msg << "Alquimia_Chemistry_PK::XMLParameters(): \n";
    msg << "  No 'geochemical conditions' list was found in 'State->initial conditions'!\n";
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::ParameterList geochem_conditions = initial_conditions.sublist("geochemical conditions");
  ParseChemicalConditionRegions(geochem_conditions, chem_initial_conditions_);
  if (chem_initial_conditions_.empty())
  {
    msg << "Alquimia_Chemistry_PK::XMLParameters(): \n";
    msg << "  No geochemical conditions were found in 'State->initial conditions->geochemical conditions'!\n";
    Exceptions::amanzi_throw(msg);
  }

  // Other settings.
  set_max_time_step(chem_param_list_.get<double>("Max Time Step (s)", 9.9e9));

}  // end XMLParameters()

void Alquimia_Chemistry_PK::CopyAmanziStateToAlquimia(const int cell_id,
                                                      Teuchos::RCP<const Epetra_MultiVector> aqueous_components,
                                                      AlquimiaMaterialProperties& mat_props,
                                                      AlquimiaState& state,
                                                      AlquimiaAuxiliaryData& aux_data)
{
  chemistry_state_->CopyToAlquimia(cell_id, aqueous_components, mat_props, state, aux_data);
}  // end CopyAmanziStateToAlquimia()

void Alquimia_Chemistry_PK::CopyAlquimiaStateToAmanzi(const int cell_id,
                                                      const AlquimiaMaterialProperties& mat_props,
                                                      const AlquimiaState& state,
                                                      const AlquimiaAuxiliaryData& aux_data,
                                                      const AlquimiaAuxiliaryOutputData& aux_output,
                                                      Teuchos::RCP<const Epetra_MultiVector> aqueous_components)
{
  chemistry_state_->CopyFromAlquimia(cell_id, mat_props, state, aux_data, aqueous_components);

  // Auxiliary output.
  if (aux_output_ != Teuchos::null) 
  {
    std::vector<std::string> mineralNames, primaryNames;
    chem_engine_->GetMineralNames(mineralNames);
    chem_engine_->GetPrimarySpeciesNames(primaryNames);
    int numAqueousComplexes = chem_engine_->NumAqueousComplexes();
    for (unsigned int i = 0; i < aux_names_.size(); i++) 
    {
      if (aux_names_.at(i) == "pH") 
      {
        double* cell_aux_output = (*aux_output_)[i];
        cell_aux_output[cell_id] = aux_output.pH;
      }
      else if (aux_names_.at(i).find("mineral_saturation_index") != std::string::npos)
      {
        for (int j = 0; j < mineralNames.size(); ++j)
        {
          std::string full_name = string("mineral_saturation_index_") + mineralNames[j];
          if (aux_names_.at(i) == full_name)
          {
            double* cell_aux_output = (*aux_output_)[i];
            cell_aux_output[cell_id] = aux_output.mineral_saturation_index.data[j];
          }
        }
      }
      else if (aux_names_.at(i).find("mineral_reaction_rate") != std::string::npos)
      {
        for (int j = 0; j < mineralNames.size(); ++j)
        {
          std::string full_name = string("mineral_reaction_rate_") + mineralNames[j];
          if (aux_names_.at(i) == full_name)
          {
            double* cell_aux_output = (*aux_output_)[i];
            cell_aux_output[cell_id] = aux_output.mineral_reaction_rate.data[j];
          }
        }
      }
      else if (aux_names_.at(i) == "primary_free_ion_concentration")
      {
        for (int j = 0; j < primaryNames.size(); ++j)
        {
          std::string full_name = string("primary_free_ion_concentration_") + primaryNames[j];
          if (aux_names_.at(i) == full_name)
          {
            double* cell_aux_output = (*aux_output_)[i];
            cell_aux_output[cell_id] = aux_output.primary_free_ion_concentration.data[j];
          }
        }
      }
      else if (aux_names_.at(i) == "primary_activity_coeff")
      {
        for (int j = 0; j < primaryNames.size(); ++j)
        {
          std::string full_name = string("primary_activity_coeff_") + primaryNames[j];
          if (aux_names_.at(i) == full_name)
          {
            double* cell_aux_output = (*aux_output_)[i];
            cell_aux_output[cell_id] = aux_output.primary_activity_coeff.data[j];
          }
        }
      }
      else if (aux_names_.at(i) == "secondary_free_ion_concentration")
      {
        for (int j = 0; j < numAqueousComplexes; ++j)
        {
          char num_str[16];
          snprintf(num_str, 15, "%d", j);
          std::string full_name = string("secondary_free_ion_concentration_") + string(num_str);
          if (aux_names_.at(i) == full_name)
          {
            double* cell_aux_output = (*aux_output_)[i];
            cell_aux_output[cell_id] = aux_output.secondary_free_ion_concentration.data[j];
          }
        }
      }
      else if (aux_names_.at(i) == "secondary_activity_coeff")
      {
        for (int j = 0; j < numAqueousComplexes; ++j)
        {
          char num_str[16];
          snprintf(num_str, 15, "%d", j);
          std::string full_name = string("secondary_activity_coeff_") + string(num_str);
          if (aux_names_.at(i) == full_name)
          {
            double* cell_aux_output = (*aux_output_)[i];
            cell_aux_output[cell_id] = aux_output.secondary_activity_coeff.data[j];
          }
        }
      }
    }
  }
}  // end CopyAlquimiaStateToAmanzi()

/*******************************************************************************
 **
 **  MPC interface functions
 **
 ******************************************************************************/
Teuchos::RCP<Epetra_MultiVector> Alquimia_Chemistry_PK::get_total_component_concentration(void) const 
{
  return chemistry_state_->total_component_concentration();
}  // end get_total_component_concentration()


// This helper advances the solution on a single cell within Amanzi's state.
// It returns the number of iterations taken to obtain the advanced solution, 
// or -1 if an error occurred.
int Alquimia_Chemistry_PK::AdvanceSingleCell(double delta_time, 
                                             Teuchos::RCP<const Epetra_MultiVector> total_component_concentration_star,
                                             int cellIndex)
{
  // Copy the state and material information from Amanzi's state within 
  // this cell to Alquimia.
  CopyAmanziStateToAlquimia(cellIndex, 
                            total_component_concentration_star, 
                            alq_mat_props_, alq_state_, alq_aux_data_);

  // Do the reaction.
  int num_iterations;
  chem_engine_->Advance(delta_time, alq_mat_props_, alq_state_, 
                        alq_aux_data_, alq_aux_output_, num_iterations);

  // Move the information back into Amanzi's state, updating the given total concentration vector.
  CopyAlquimiaStateToAmanzi(cellIndex, alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_, 
                            chemistry_state_->total_component_concentration());

  return num_iterations;
}

/*******************************************************************************
 **
 ** Alquimia_Chemistry_PK::advance()
 **
 ** Notes:
 **
 **   - the MPC will call this function to advance the state with this
 ** particular process kernel
 **
 **   - this is how to get the total component concentration
 ** CS->get_total_component_concentration()
 **
 **   - please update the argument to this function called tcc_star
 ** with the result of your chemistry computation which is the total
 ** component concentration ^star
 **
 **   - see the Chemistry_State for the other available data in the
 ** chemistry specific state
 **
 *******************************************************************************/

void Alquimia_Chemistry_PK::advance(
    const double& delta_time,
    Teuchos::RCP<const Epetra_MultiVector> total_component_concentration_star) 
{
  Errors::Message msg;

  current_time_ = saved_time_ + delta_time;

  // shorter name for the state that came out of transport
  Teuchos::RCP<const Epetra_MultiVector> tcc_star =
      total_component_concentration_star;


  // TODO(bandre): use size of the porosity vector as indicator of size for
  // now... should get data from the mesh...?
  int num_cells = chemistry_state_->porosity()->MyLength();

  int max_iterations = 0;
  int imax = -999;
  int ave_iterations = 0;
  int imin = -999;
  int min_iterations = 10000000;

  // Now loop through all the regions and advance the chemistry.
  int ierr = 0;
  for (int cell = 0; cell < num_cells; ++cell)
  {
    int num_iterations = AdvanceSingleCell(delta_time, tcc_star, cell);
    if (num_iterations > 0)
    {
      //std::stringstream message;
      //message << "--- " << cell << "\n";
      //beaker_components_.Display(message.str().c_str());
      if (max_iterations < num_iterations) 
      {
        max_iterations = num_iterations;
        imax = cell;
      }
      if (min_iterations > num_iterations) 
      {
        min_iterations = num_iterations;
        imin = cell;
      }
      ave_iterations += num_iterations;
    } 
    else
    {
      ierr = 1;
    }
  }

}  // end advance()


// the MPC will call this function to signal to the
// process kernel that it has accepted the
// state update, thus, the PK should update
// possible auxilary state variables here
void Alquimia_Chemistry_PK::commit_state(Teuchos::RCP<Chemistry_State> chem_state,
                                         const double& delta_time) 
{
  saved_time_ += delta_time;
}  // end commit_state()



Teuchos::RCP<Epetra_MultiVector> Alquimia_Chemistry_PK::get_extra_chemistry_output_data() 
{
  // This vector is updated during the initialization and advance of 
  // the geochemistry, so we simply return it here.
  return aux_output_;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
