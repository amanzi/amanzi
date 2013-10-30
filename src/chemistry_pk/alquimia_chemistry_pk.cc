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

// For now, we use the "interim" chemistry input spec, in which the 
// initial and boundary conditions are specified within the Chemistry
// parameter list.
#define USE_INTERIM_INPUT_SPEC 1

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
 **
 **    where foo refers to a vector of component concentrations for a
 **    single component.
 **
 **
 ******************************************************************************/

Alquimia_Chemistry_PK::Alquimia_Chemistry_PK(const Teuchos::ParameterList& param_list,
                                             Teuchos::RCP<Chemistry_State> chem_state,
                                             Teuchos::RCP<Chemistry_Engine> chemistry_engine)
    : debug_(false),
      display_free_columns_(false),
      max_time_step_(9.9e9),
      chemistry_state_(chem_state),
      main_param_list_(param_list),
      chem_param_list_(),
      chem_initialized_(false),
      chem_engine_(chemistry_engine),
      current_time_(0.0),
      saved_time_(0.0) 
{
  // We need the top-level parameter list.
//  assert(param_list.name() == "Main");
  chem_param_list_ = main_param_list_.sublist("Chemistry");
}  // end Alquimia_Chemistry_PK()

Alquimia_Chemistry_PK::~Alquimia_Chemistry_PK() 
{
}  // end ~Alquimia_Chemistry_PK()

// This helper performs initialization on a single cell within Amanzi's state.
// It returns an error code that indicates success (0) or failure (1).
int Alquimia_Chemistry_PK::InitializeSingleCell(int cellIndex, const std::string& condition) 
{
  // Here, we create arrays for transferring state information to the chemistry engine.
  int N_primary = chem_engine_->NumPrimarySpecies();
  int N_sorbed = chem_engine_->NumSorbedSpecies();
  int N_minerals = chem_engine_->NumMinerals();
  int N_surf_sites = chem_engine_->NumSurfaceSites();
  int N_ion_exch = chem_engine_->NumIonExchangeSites();
  int N_free_ions = chem_engine_->NumFreeIonSpecies();
  int N_isotherm_species = chem_engine_->NumIsothermSpecies();
  std::vector<double> component_concentrations(N_primary), sorbed_components(N_sorbed),
                      mineral_volume_fractions(N_minerals), mineral_specific_surface_areas(N_minerals),
                      cation_exchange_capacity(N_ion_exch), sorption_sites(N_surf_sites),
                      free_ion_species(N_free_ions), primary_activity_coeffs(N_primary), // FIXME: correct size?
                      secondary_activity_coeffs(N_primary), // FIXME: correct size?
                      ion_exchange_ref_cation_concs(N_ion_exch), // FIXME: correct size?
                      surface_complex_free_site_concs(N_surf_sites), // FIXME: correct size?
                      isotherm_kd(N_isotherm_species), 
                      isotherm_freundlich_n(N_isotherm_species),
                      isotherm_langmuir_b(N_isotherm_species);
  double water_density, porosity, volume, saturation, pH;

  // Copy the state and material information from Amanzi's state within 
  // this cell to Alquimia.
  CopyAmanziStateToBuffers(cellIndex, 
                           chemistry_state_->total_component_concentration(), 
                           component_concentrations,
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

  // Do the initializætion.
  chem_engine_->EnforceCondition(condition,
                                 current_time_,
                                 &component_concentrations[0],
                                 &sorbed_components[0],
                                 &mineral_volume_fractions[0],
                                 &mineral_specific_surface_areas[0],
                                 &cation_exchange_capacity[0],
                                 &sorption_sites[0],
                                 water_density,
                                 porosity,
                                 volume,
                                 saturation,
                                 &isotherm_kd[0],
                                 &isotherm_freundlich_n[0],
                                 &isotherm_langmuir_b[0],
                                 &free_ion_species[0],
                                 &primary_activity_coeffs[0],
                                 &secondary_activity_coeffs[0],
                                 &ion_exchange_ref_cation_concs[0],
                                 &surface_complex_free_site_concs[0],
                                 pH);

  // Move the information back into Amanzi's state.
  CopyBuffersToAmanziState(cellIndex, 
                           component_concentrations,
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

  return 0;
}

void Alquimia_Chemistry_PK::InitializeChemistry(void) 
{
  Errors::Message msg;
  if (debug()) 
  {
    std::cout << "  Alquimia_Chemistry_PK::InitializeChemistry()" << std::endl;
  }

  // Read XML parameters from our input file.
  XMLParameters();

  chemistry_state_->AllocateAdditionalChemistryStorage(number_aqueous_components());

  // Make storage within Amanzi for auxiliary output.
  SetupAuxiliaryOutput();

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
      if (debug()) 
      {
        std::cout << "Initial speciation in cell " << cell << std::endl;
      }
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

//  chem_out->Write(kVerbose, "ChemistryPK::InitializeChemistry(): initialization was successful.\n");

  chem_initialized_ = true;
}  // end InitializeChemistry()

/*******************************************************************************
 **
 **  initialization helper functions
 **
 ******************************************************************************/
void Alquimia_Chemistry_PK::ParseChemicalConditions(const Teuchos::ParameterList& param_list,
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
 
    // At the moment, we only support constraints supplied through the backend engine.
    std::string engine_constraint = chem_engine_->Name() + std::string(" Constraint");
    if (!cond_sublist.isType<std::string>(engine_constraint))
    {
      msg << "Chemistry_PK::XMLParameters(): \n";
      msg << "  Geochemical condition '" << cond_name << "' has no '" << engine_constraint << "' entry.\n";
      msg << "  (Must be of type Array(string).)\n";
      Exceptions::amanzi_throw(msg);
    }
    std::string condition = cond_sublist.get<std::string>(engine_constraint);

    // Apply this condition to all desired regions.
    if (!cond_sublist.isType<Teuchos::Array<std::string> >("Assigned Regions"))
    {
      msg << "Chemistry_PK::XMLParameters(): \n";
      msg << "  Geochemical condition '" << cond_name << "' has no valid 'Assigned Regions' entry.\n";
      Exceptions::amanzi_throw(msg);
    }
    Teuchos::Array<std::string> regions = cond_sublist.get<Teuchos::Array<std::string> >("Assigned Regions");
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = chemistry_state_->mesh_maps();
    for (size_t r = 0; r < regions.size(); ++r)
    {
      // We allow for cell-based and face-based regions to accommodate both 
      // initial and boundary conditions.
      if (!mesh->valid_set_name(regions[r], AmanziMesh::CELL) &&
          !mesh->valid_set_name(regions[r], AmanziMesh::FACE))
      {
        msg << "Chemistry_PK::XMLParameters(): \n";
        msg << "  Invalid region '" << regions[r] << "' given for geochemical condition '" << cond_name << "'.\n";
        Exceptions::amanzi_throw(msg);
      }
      conditions[regions[r]] = condition;
    }
  }
}

void Alquimia_Chemistry_PK::XMLParameters(void) 
{
  Errors::Message msg;

  // NOTE that our parameter list should be the top-level parameter list "Main", 
  // not the "Chemistry" one used by the native Amanzi Chemistry PK.

  // Debugging flag.
  if (chem_param_list_.isParameter("Debug Process Kernel")) 
  {
    set_debug(true);
  }

  // Verbosity.
//  if (chem_param_list_.isParameter("Verbosity")) 
//  {
//    Teuchos::Array<std::string> verbosity_list = chem_param_list_.get<Teuchos::Array<std::string> >("Verbosity");
//    Teuchos::Array<std::string>::const_iterator name;
//    for (name = verbosity_list.begin(); name != verbosity_list.end(); ++name) 
//    {
//      chem_out->AddLevel(*name);
//    }
//  }

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
//        message << "ChemistryPK::XMLParameters(): unknown value in 'Auxiliary Data' list: " 
//                << *name << std::endl;
//        chem_out->Write(kWarning, message);
//      }
    }
  }

  // Now set up chemical conditions based on initial and boundary conditions
  // elsewhere in the file.

  // Initial conditions.
  if (!main_param_list_.isSublist("State"))
//  if (!main_param_list_.isSublist("Initial Conditions"))
  {
    msg << "Chemistry_PK::XMLParameters(): \n";
    msg << "  No 'State' sublist was found!\n";
    Exceptions::amanzi_throw(msg);
  }
#if USE_INTERIM_INPUT_SPEC
  Teuchos::ParameterList initial_conditions = chem_param_list_.sublist("Initial Conditions");
#else
  Teuchos::ParameterList state_list = main_param_list_.sublist("state");
  Teuchos::ParameterList initial_conditions = state_list.sublist("initial conditions");
#endif
  ParseChemicalConditions(initial_conditions, chem_initial_conditions_);
  if (chem_initial_conditions_.empty())
  {
    msg << "Chemistry_PK::XMLParameters(): \n";
    msg << "  No chemical conditions were found in 'initial conditions'!\n";
    Exceptions::amanzi_throw(msg);
  }

  // Other settings.
  set_max_time_step(chem_param_list_.get<double>("Max Time Step (s)", 9.9e9));

}  // end XMLParameters()

void Alquimia_Chemistry_PK::SetupAuxiliaryOutput(void) 
{
  // requires that Alquimia's Setup method has already been called!
  if (debug()) 
  {
    std::cout << "  Alquimia_Chemistry_PK::SetupAuxiliaryOutput()" << std::endl;
  }
  // TODO(bandre): this indexing scheme will not be appropriate when
  // additional types of aux data are requested, e.g. mineral SI.....
  unsigned int nvars = aux_names_.size();
  std::string name;
  for (unsigned int i = 0; i < nvars; i++) 
  {
    if (aux_names_.at(i) == "pH") 
    {
      name = "H+";
    } 
    else 
    {
      Errors::Message msg;
      msg << "ChemistryPK::SetupAuxiliaryOutput() : "
          << "Output was requested for '" << aux_names_.at(i) 
          << "' (" << name
          << ") but no auxiliary chemistry varibles of this name were found.\n";
      Exceptions::amanzi_throw(msg);
    } 
  }

  // create the Epetra_MultiVector that will hold the data
  if (nvars > 0) 
  {
    aux_data_ =
        Teuchos::rcp(new Epetra_MultiVector(
            chemistry_state_->mesh_maps()->cell_map(false), nvars));
  } 
  else 
  {
    aux_data_ = Teuchos::null;
  }
}  // end SetupAuxiliaryOutput()

void Alquimia_Chemistry_PK::CopyAmanziStateToBuffers(const int cell_id,
                                                     Teuchos::RCP<const Epetra_MultiVector> aqueous_components, 
                                                     std::vector<double>& component_concentrations,
                                                     std::vector<double>& sorbed_components,
                                                     std::vector<double>& mineral_volume_fractions,
                                                     std::vector<double>& mineral_specific_surface_areas,
                                                     std::vector<double>& cation_exchange_capacity,
                                                     std::vector<double>& sorption_sites,
                                                     double& water_density,
                                                     double& porosity,
                                                     double& volume,
                                                     double& saturation,
                                                     std::vector<double>& isotherm_kd,
                                                     std::vector<double>& isotherm_freundlich_n,
                                                     std::vector<double>& isotherm_langmuir_b)
{
  water_density = (*chemistry_state_->water_density())[cell_id];
  porosity = (*chemistry_state_->porosity())[cell_id];

  for (unsigned int c = 0; c < number_aqueous_components(); c++) 
  {
    double* cell_components = (*aqueous_components)[c];
    double total_mobile = cell_components[cell_id];
    component_concentrations[c] = total_mobile;
    if (using_sorption()) 
    {
      double* cell_total_sorbed = (*chemistry_state_->total_sorbed())[c];
      double immobile = cell_total_sorbed[cell_id];
      sorbed_components[c] = immobile;
    }  // end if(using_sorption)
  }

  // minerals
  for (unsigned int m = 0; m < number_minerals(); m++) 
  {
    double* cell_minerals = (*chemistry_state_->mineral_volume_fractions())[m];
    mineral_volume_fractions[m] = cell_minerals[cell_id];
    if (chemistry_state_->mineral_specific_surface_area() != Teuchos::null) 
    {
      double* cells_ssa = (*chemistry_state_->mineral_specific_surface_area())[m];
      mineral_specific_surface_areas[m] = cells_ssa[cell_id];
    }
  }

  // ion exchange
  if (number_ion_exchange_sites() > 0) 
  {
    for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) 
    {
      double* cell_ion_exchange_sites = (*chemistry_state_->ion_exchange_sites())[i];
      cation_exchange_capacity[i] = cell_ion_exchange_sites[cell_id];
    }
  }
  
  // surface complexation
  if (number_sorption_sites() > 0) 
  {
    for (int s = 0; s < number_sorption_sites(); ++s) 
    {
      // FIXME: Need site density names, too?
      double* cell_sorption_sites = (*chemistry_state_->sorption_sites())[s];
      sorption_sites[s] = cell_sorption_sites[cell_id];
      // TODO(bandre): need to save surface complexation free site conc here!
    }
  }

  // Material properties.
  volume = (*chemistry_state_->volume())[cell_id];
  saturation = (*chemistry_state_->water_saturation())[cell_id];

  // sorption isotherms
  if (using_sorption_isotherms()) 
  {
    for (unsigned int i = 0; i < number_aqueous_components(); ++i) 
    {
      double* cell_data = (*chemistry_state_->isotherm_kd())[i];
      isotherm_kd[i] = cell_data[cell_id];
      
      cell_data = (*chemistry_state_->isotherm_freundlich_n())[i];
      isotherm_freundlich_n[i] = cell_data[cell_id];
      
      cell_data = (*chemistry_state_->isotherm_langmuir_b())[i];
      isotherm_langmuir_b[i] = cell_data[cell_id];
    }
  }

}  // end CopyAmanziStateToBuffers()

void Alquimia_Chemistry_PK::CopyBuffersToAmanziState(const int cell_id,
                                                     const std::vector<double>& component_concentrations,
                                                     const std::vector<double>& sorbed_components,
                                                     const std::vector<double>& mineral_volume_fractions,
                                                     const std::vector<double>& mineral_specific_surface_areas,
                                                     const std::vector<double>& cation_exchange_capacity,
                                                     const std::vector<double>& sorption_sites,
                                                     double water_density,
                                                     double porosity,
                                                     double volume,
                                                     double saturation,
                                                     const std::vector<double>& isotherm_kd,
                                                     const std::vector<double>& isotherm_freundlich_n,
                                                     const std::vector<double>& isotherm_langmuir_b,
                                                     const std::vector<double>& free_ion_species,
                                                     const std::vector<double>& primary_activity_coeffs,
                                                     const std::vector<double>& secondary_activity_coeffs,
                                                     const std::vector<double>& ion_exchange_ref_cation_concs,
                                                     const std::vector<double>& surface_complex_free_site_concs,
                                                     double pH)
{
  // If the chemistry has modified the porosity and/or density, it needs to 
  // be updated here.
  // NOTE: At the moment, these accessors are const-only in Chemistry_State.
  //(*chemistry_state_->water_density())[cell_id] = water_density;
  //(*chemistry_state_->porosity())[cell_id] = porosity;

  for (unsigned int c = 0; c < number_aqueous_components(); c++) 
  {
    double mobile = component_concentrations[c];

    // NOTE: We place our mobile concentrations into the chemistry state's
    // NOTE: total component concentrations, not the aqueous concentrations.
    // NOTE: I assume this has to do with our operator splitting technique.
    double total = mobile;
    double* cell_components = (*chemistry_state_->total_component_concentration())[c];
    cell_components[cell_id] = total;

    if (using_sorption())
    {
      double immobile = sorbed_components[c];
      double* cell_total_sorbed = (*chemistry_state_->total_sorbed())[c];
      cell_total_sorbed[cell_id] = immobile;
    }
  }

  // minerals
  for (unsigned int m = 0; m < number_minerals(); m++) 
  {
    double* cell_minerals = (*chemistry_state_->mineral_volume_fractions())[m];
    cell_minerals[cell_id] = mineral_volume_fractions[m];
    if (chemistry_state_->mineral_specific_surface_area() != Teuchos::null) 
    {
      cell_minerals = (*chemistry_state_->mineral_specific_surface_area())[m];
      cell_minerals[cell_id] = mineral_specific_surface_areas[m];
    }
  }

  // ion exchange
  for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) 
  {
    double* cell_ion_exchange_sites = (*chemistry_state_->ion_exchange_sites())[i];
    cell_ion_exchange_sites[cell_id] = cation_exchange_capacity[i];
  }

  // surface complexation
  if (number_sorption_sites() > 0)
  {
    for (unsigned int i = 0; i < number_sorption_sites(); i++) 
    {
      double* cell_sorption_sites = (*chemistry_state_->sorption_sites())[i];
      cell_sorption_sites[cell_id] = sorption_sites[i];
    }
  }

  // free ion concentrations
  int offset = 0;
  for (int i = 0; i < number_aqueous_components(); ++i, ++offset) 
  {
    double* cell_free_ion = (*chemistry_state_->free_ion_species())[i];
    cell_free_ion[cell_id] = free_ion_species[i];
  }

  // activity coefficients
  for (int i = 0; i < number_aqueous_components(); ++i, ++offset) 
  {
    double* cells = (*chemistry_state_->primary_activity_coeff())[i];
    cells[cell_id] = primary_activity_coeffs[i];
  }
  for (int i = 0; i < number_aqueous_components(); ++i, ++offset) 
  {
    double* cells = (*chemistry_state_->secondary_activity_coeff())[i];
    cells[cell_id] = secondary_activity_coeffs[i];
  }

  // ion exchange ref cation concentrations
  for (int i = 0; i < chemistry_state_->number_of_ion_exchange_sites(); ++i, ++offset) 
  {
    double* cells = (*chemistry_state_->ion_exchange_ref_cation_conc())[i];
    cells[cell_id] = ion_exchange_ref_cation_concs[i];
  }

  // surface complexation
  if (using_sorption())
  {
    for (int i = 0; i < chemistry_state_->number_of_sorption_sites(); ++i, ++offset) 
    {
      double* cells = (*chemistry_state_->surface_complex_free_site_conc())[i];
      cells[cell_id] = surface_complex_free_site_concs[i];
    }
  }

  // Auxiliary output.
  if (aux_data_ != Teuchos::null) 
  {
    for (unsigned int i = 0; i < aux_names_.size(); i++) 
    {
      if (aux_names_.at(i) == "pH") 
      {
        double* cell_aux_data = (*aux_data_)[i];
        cell_aux_data[cell_id] = pH;
      }
      // FIXME: For now, we don't support mineral saturation index or 
      // FIXME: mineral reaction rate, but when we do, these need to 
      // FIXME: be updated here as well.
    }
  }

  if (using_sorption_isotherms()) 
  {
    for (unsigned int i = 0; i < number_aqueous_components(); ++i) 
    {
      double* cell_data = (*chemistry_state_->isotherm_kd())[i];
      cell_data[cell_id] = isotherm_kd[i];

      cell_data = (*chemistry_state_->isotherm_freundlich_n())[i];
      cell_data[cell_id] = isotherm_freundlich_n[i];

      cell_data = (*chemistry_state_->isotherm_langmuir_b())[i];
      cell_data[cell_id] = isotherm_langmuir_b[i];
    }
  }

}  // end CopyBuffersToAmanziState()

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
  // Here, we create arrays for transferring state information to the chemistry engine.
  int N_primary = chem_engine_->NumPrimarySpecies();
  int N_sorbed = chem_engine_->NumSorbedSpecies();
  int N_minerals = chem_engine_->NumMinerals();
  int N_surf_sites = chem_engine_->NumSurfaceSites();
  int N_ion_exch = chem_engine_->NumIonExchangeSites();
  int N_free_ions = chem_engine_->NumFreeIonSpecies();
  int N_isotherm_species = chem_engine_->NumIsothermSpecies();
  std::vector<double> component_concentrations(N_primary), sorbed_components(N_sorbed),
                      mineral_volume_fractions(N_minerals), mineral_specific_surface_areas(N_minerals),
                      cation_exchange_capacity(N_ion_exch), sorption_sites(N_surf_sites),
                      free_ion_species(N_free_ions), primary_activity_coeffs(N_primary), // FIXME: correct size?
                      secondary_activity_coeffs(N_primary), // FIXME: correct size?
                      ion_exchange_ref_cation_concs(N_ion_exch), // FIXME: correct size?
                      surface_complex_free_site_concs(N_surf_sites), // FIXME: correct size?
                      isotherm_kd(N_isotherm_species), 
                      isotherm_freundlich_n(N_isotherm_species),
                      isotherm_langmuir_b(N_isotherm_species);
  double water_density, porosity, volume, saturation, pH;

  // Copy the state and material information from Amanzi's state within 
  // this cell to Alquimia.
  CopyAmanziStateToBuffers(cellIndex, 
                           total_component_concentration_star, 
                           component_concentrations,
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

  // Do the reaction.
  int num_iterations;
  chem_engine_->Advance(delta_time,
                        &component_concentrations[0],
                        &sorbed_components[0],
                        &mineral_volume_fractions[0],
                        &mineral_specific_surface_areas[0],
                        &cation_exchange_capacity[0],
                        &sorption_sites[0],
                        water_density,
                        porosity,
                        volume,
                        saturation,
                        &isotherm_kd[0],
                        &isotherm_freundlich_n[0],
                        &isotherm_langmuir_b[0],
                        &free_ion_species[0],
                        &primary_activity_coeffs[0],
                        &secondary_activity_coeffs[0],
                        &ion_exchange_ref_cation_concs[0],
                        &surface_complex_free_site_concs[0],
                        pH,
                        num_iterations);

  // Move the information back into Amanzi's state.
  CopyBuffersToAmanziState(cellIndex, 
                           component_concentrations,
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

  // This shouldn't be necessary, as the Chemistry Engine handles errors.
#if 0
  int recv = 0;
  chemistry_state_->mesh_maps()->get_comm()->MaxAll(&ierr, &recv, 1);
  if (recv != 0) 
  {
    msg << "Error in Alquimia_Chemistry_PK::advance: ";
    msg << chem_status_.message;
    Exceptions::amanzi_throw(msg); 
  }
#endif
  
#ifdef GLENN_DEBUG
  if (debug() == kDebugChemistryProcessKernel) 
  {
    std::stringstream message;
    message << "  Alquimia_Chemistry_PK::advance() : "
            << "max iterations - " << max_iterations << " " << "  cell id: " << imax << std::endl;
    message << "  Alquimia_Chemistry_PK::advance() : "
            << "min iterations - " << min_iterations << " " << "  cell id: " << imin << std::endl;
    message << "  Alquimia_Chemistry_PK::advance() : "
            << "ave iterations - " << static_cast<float>(ave_iterations) / num_cells << std::endl;
  }
#endif

}  // end advance()


// the MPC will call this function to signal to the
// process kernel that it has accepted the
// state update, thus, the PK should update
// possible auxilary state variables here
void Alquimia_Chemistry_PK::commit_state(Teuchos::RCP<Chemistry_State> chem_state,
                                         const double& delta_time) 
{
//  if (debug() == kDebugChemistryProcessKernel) 
//  {
//    chem_out->Write(kVerbose,
//                    "  Alquimia_Chemistry_PK::commit_state() : Committing internal state.\n");
//  }

  saved_time_ += delta_time;

  // FIXME: I'm not sure there's anything else we have to do here.
#if 0
  if (debug() && false) 
  {
    chem_->Speciate(&beaker_components_, beaker_parameters_);
    chem_->DisplayResults();
    chem_->DisplayTotalColumnHeaders(display_free_columns_);
    chem_->DisplayTotalColumns(saved_time_, beaker_components_, true);
  }
#endif
}  // end commit_state()



Teuchos::RCP<Epetra_MultiVector> Alquimia_Chemistry_PK::get_extra_chemistry_output_data() 
{
  // This vector is updated during the initialization and advance of 
  // the geochemistry, so we simply return it here.
  return aux_data_;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
