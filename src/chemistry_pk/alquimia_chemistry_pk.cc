/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "alquimia_chemistry_pk.hh"

#include <string>
#include <algorithm>

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

#include "simple_thermo_database.hh"

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
 **
 **    where foo refers to a vector of component concentrations for a
 **    single component.
 **
 **
 ******************************************************************************/

// global ChemistryOutput object in the Amanzi::chemisry
// namespace that will be used by the chemistry library
extern ChemistryOutput* chem_out;

Alquimia_Chemistry_PK::Alquimia_Chemistry_PK(const Teuchos::ParameterList& param_list,
                                             Teuchos::RCP<Chemistry_State> chem_state)
    : debug_(false),
      display_free_columns_(false),
      max_time_step_(9.9e9),
      current_time_(0.0),
      chem_initialized_(false),
      saved_time_(0.0) {

}  // end Alquimia_Chemistry_PK()

Alquimia_Chemistry_PK::~Alquimia_Chemistry_PK() {
  chem_.Shutdown(&chem_data_.engine_state, &chem_status_);
  FreeAlquimiaData(&chem_data_);
  FreeAlquimiaEngineStatus(&chem_status_);
}  // end ~Alquimia_Chemistry_PK()

// This helper performs initialization on a single cell within Amanzi's state.
// It returns an error code that indicates success (0) or failure (1).
int Alquimia_Chemistry_PK::InitializeSingleCell(int cellIndex) {

  // Copy the state and material information from Amanzi's state within 
  // this cell to Alquimia.
  CopyAmanziStateToAlquimia(cell, chemistry_state_->total_component_concentration());
  CopyAmanziMaterialPropertiesToAlquimia(cell, chemistry_state_->total_component_concentration());

  // Apply the geochemical conditions.
  int ierr = 0;
  try {
    chem.ProcessCondition(&chem_data_.engine_state,
                          &(alquimia_conditions_.data[i]),
                          &chem_data_.material_properties,
                          &chem_data_.state,
                          &chem_data_.aux_data,
                          &chem_status_);
      //chem_->DisplayTotalColumns(static_cast<double>(-cell), beaker_components_, display_free_columns_);
      // solve for initial free-ion concentrations
//      chem_->Speciate(&beaker_components_, beaker_parameters_);
//      CopyBeakerStructuresToCellState(cell);
  } catch (AmanziException& geochem_error) {
    ierr = 1;
    // std::cout << "ChemistryPK::InitializeChemistry(): " 
    //           << "An error occured while initializing chemistry in cell: " 
    //           << cell << ".\n" << geochem_error.what() << std::endl;
    // beaker_components_.Display("-- input components:\n");
    // Exceptions::amanzi_throw(geochem_error);
  }

  // Retrieve the auxiliary output data.
  chem_.GetAuxiliaryOutput(&chem_data_.engine_state,
                           &chem_data_.material_properties,
                           &chem_data_.state,
                           &chem_data_.aux_data,
                           &chem_data_.aux_output,
                           &chem_status_);

  // Copy the state information back out.
  CopyAlquimiaStateToAmanzi(cell);
  CopyAlquimiaMaterialPropertiesToAmanzi(cell);

  return ierr;
}

void Alquimia_Chemistry_PK::InitializeChemistry(void) {
  if (debug()) {
    std::cout << "  Alquimia_Chemistry_PK::InitializeChemistry()" << std::endl;
  }

  // Read XML parameters from our input file.
  // FIXME: Not implemented properly. Among other things, needs to set the 
  // FIXME: input file for Alquimia's "engine."
  XMLParameters();

  // We can't initialize chemistry twice in a row.
  if (chem_initialized_)
  {
    AmanziException geochem_error("Alquimia_Chemistry_PK: initialized chemistry twice!");
    Exceptions::amanzi_throw(geochem_error); 
  }

  // All alquimia function calls require a status object.
  AllocateAlquimiaEngineStatus(&chem_status_);
  CreateAlquimiaInterface(chem_engine_name_.c_str(), &chem_, &chem_status_);
  if (chem_status_.error != 0) {
    std::cout << chem_status_.message << std::endl;
    AmanziException geochem_error("Alquimia_Chemistry_PK: Could not create an interface to Alquimia.");
    Exceptions::amanzi_throw(geochem_error); 
  }

  // Set up Alquimia, get sizes for data.
  chem.Setup(chem_engine_inputfile_.c_str(),
             &chem_data_.engine_state,
             &chem_data_.sizes,
             &chem_data_.functionality,
             &chem_status_);
  if (chem_status_.error != 0) {
    std::cout << chem_status_.message << std::endl;
    PrintAlquimiaSizes(&chem_data_.sizes);
    AmanziException geochem_error("Error in Alquimia_Chemistry_PK::InitializeChemistry 1");
    Exceptions::amanzi_throw(geochem_error); 
  }

  // Allocate data storage now that we have our space requirements.
  AllocateAlquimiaData(&chem_data_);

  // Get the problem metadata (species and mineral names, etc).
  chem_.GetProblemMetaData(&chem_data_.engine_state,
                           &chem_data_.meta_data,
                           &chem_status_);
  if (chem_status_.error != 0) {
    std::cout << chem_status_.message << std::endl;
    PrintAlquimiaProblemMetaData(&chem_data_.meta_data);
    AmanziException geochem_error("Error in Alquimia_Chemistry_PK::InitializeChemistry 1");
    Exceptions::amanzi_throw(geochem_error); 
  }
  // FIXME: Verify that sizes/ordering are consistent with MPC state.

  // Initialize Alquimia's state and material properties with
  // appropriate values from the our cell properties.

  // NOTE: we need to perform the initialization on the first cell to determine 
  // NOTE: the number of secondary activity coefficients, which we'll stash in 
  // NOTE: the chemistry state.
  int ierr = InitializeSingleCell(0);
  if (ierr != 0)
  {
    AmanziException geochem_error("Error in Alquimia_Chemistry_PK::InitializeChemistry 1");
    Exceptions::amanzi_throw(geochem_error); 
  }
  chemistry_state_->AllocateAdditionalChemistryStorage(number_aqueous_components());

  SetupAuxiliaryOutput();

  // Now loop through all the cells and initialize
  chem_out->Write(kVerbose, "ChemistryPK: Initializing chemistry in all cells...\n");
  int num_cells = chemistry_state_->porosity()->MyLength();
  ierr = 0;
  for (cell = 0; cell < num_cells; ++cell) {
    if (debug()) {
      std::cout << "Initial speciation in cell " << cell << std::endl;
    }

    ierr = InitializeSingleCell(cell);
  }
  recv = 0;
  // figure out if any of the processes threw an error, if so all processes will re-throw
  chemistry_state_->mesh_maps()->get_comm()->MaxAll(&ierr, &recv, 1);
  if (recv != 0) {
    AmanziException geochem_error("Error in Alquimia_Chemistry_PK::InitializeChemistry 1");
    Exceptions::amanzi_throw(geochem_error); 
  }  

  chem_out->Write(kVerbose, "ChemistryPK::InitializeChemistry(): initialization was successful.\n");

  chem_initialized_ = true;
}  // end InitializeChemistry()

/*******************************************************************************
 **
 **  initialization helper functions
 **
 ******************************************************************************/
void Alquimia_Chemistry_PK::XMLParameters(void) {
  // extract parameters from the xml list and set in the parameters
  // structure
  if (parameter_list_.isParameter("Debug Process Kernel")) {
    set_debug(true);
  }

  if (parameter_list_.isParameter("Verbosity")) {
    Teuchos::Array<std::string> verbosity_list = parameter_list_.get<Teuchos::Array<std::string> >("Verbosity");
    Teuchos::Array<std::string>::const_iterator name;
    for (name = verbosity_list.begin(); name != verbosity_list.end(); ++name) {
      chem_out->AddLevel(*name);
    }
  }
  //--------------------------------------------------------------------------
  //
  // thermo file name and format, then create the database!
  //
  //--------------------------------------------------------------------------
  if (parameter_list_.isSublist("Thermodynamic Database")) {
    Teuchos::ParameterList& tdb_list_ = parameter_list_.sublist("Thermodynamic Database");
    // get format
    // currently we only support the simple format..
    if (tdb_list_.isParameter("Format")) {
      std::string database_format = tdb_list_.get<std::string>("Format");
      if (database_format == "simple") {
        chem_ = new SimpleThermoDatabase();
      } else {
        // invalid database format...
        std::ostringstream error_stream;
        error_stream << AmanziException::kChemistryError;
        error_stream << "Chemistry_PK::XMLParameters(): \n";
        error_stream << "  In sublist 'Thermodynamic Database', the parameter 'Format' must be 'simple'.\n";
        Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));  
      }
    } else {
      // invalid database format...
      std::ostringstream error_stream;
      error_stream << AmanziException::kChemistryError;
      error_stream << "Chemistry_PK::XMLParameters(): \n";
      error_stream << "  In sublist 'Thermodynamic Database', the parameter 'Format' must be specified.\n";
      Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
    }
    beaker_parameters_ = chem_->GetDefaultParameters();
    // get file name
    if (tdb_list_.isParameter("File")) {
      beaker_parameters_.thermo_database_file = tdb_list_.get<std::string>("File");
    } else {
      std::ostringstream error_stream;
      error_stream << AmanziException::kChemistryError;
      error_stream << "Chemistry_PK::XMLParameters(): \n";
      error_stream << "  Input parameter 'File' in 'Thermodynamic Database' sublist must be specified.\n";
      Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));         
    }
  } else {
    std::ostringstream error_stream;
    error_stream << AmanziException::kChemistryError;
    error_stream << "Chemistry_PK::XMLParameters(): \n";
    error_stream << "  'Thermodynamic Database' sublist must be specified.\n";
    Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));    
  }
  //---------------------------------------------------------------------------
  //
  // activity model
  //
  //---------------------------------------------------------------------------
  beaker_parameters_.activity_model_name =
      parameter_list_.get<std::string>("Activity Model", "unit");
  // Pitzer virial coefficients database
  if (beaker_parameters_.activity_model_name == "pitzer-hwm") {
    if (parameter_list_.isParameter("Pitzer Database File")) {
      beaker_parameters_.pitzer_database =
          parameter_list_.get<std::string>("Pitzer Database File");
    } else {
      std::ostringstream error_stream;
      error_stream << AmanziException::kChemistryError;
      error_stream << "Chemistry_PK::XMLParameters(): \n";
      error_stream << "  Input parameter 'Pitzer Database File' must be specified if 'Activity Model' is 'pitzer-hwm'.\n";
      Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
      
    }
  }

  // --------------------------------------------------------------------------
  //
  // solver parameters
  //
  // --------------------------------------------------------------------------
  beaker_parameters_.tolerance =
      parameter_list_.get<double>("Tolerance", 1.0e-12);

  beaker_parameters_.max_iterations =
      static_cast<unsigned int>(parameter_list_.get<int>("Maximum Newton Iterations", 200));

  // --------------------------------------------------------------------------
  //
  // Auxiliary Data
  //
  // --------------------------------------------------------------------------
  aux_names_.clear();
  if (parameter_list_.isParameter("Auxiliary Data")) {
    Teuchos::Array<std::string> names = parameter_list_.get<Teuchos::Array<std::string> >("Auxiliary Data");
    for (Teuchos::Array<std::string>::const_iterator name = names.begin();
         name != names.end(); ++name) {
      if (*name == "pH") {
        aux_names_.push_back(*name);
      } else {
        std::stringstream message;
        message << "ChemistryPK::XMLParameters(): unknown value in 'Auxiliary Data' list: " 
                << *name << std::endl;
        chem_out->Write(kWarning, message);
      }
    }
  }
  // --------------------------------------------------------------------------
  //
  // misc other chemistry flags
  //
  // --------------------------------------------------------------------------
  set_max_time_step(parameter_list_.get<double>("Max Time Step (s)", 9.9e9));

}  // end XMLParameters()


void Alquimia_Chemistry_PK::SetupAuxiliaryOutput(void) {
  // requires that Alquimia's Setup method has already been called!
  if (debug()) {
    std::cout << "  Alquimia_Chemistry_PK::SetupAuxiliaryOutput()" << std::endl;
  }
  // TODO(bandre): this indexing scheme will not be appropriate when
  // additional types of aux data are requested, e.g. mineral SI.....
  unsigned int nvars = aux_names_.size();
  std::string name;
  aux_index_.clear();
  for (unsigned int i = 0; i < nvars; i++) {
    if (aux_names_.at(i) == "pH") {
      name = "H+";
    } else {
      name = aux_names_.at(i);
    }
    int index = chem_->GetPrimaryIndex(name);
    if (index == -1) {
        // check to make sure it is not -1, an invalid name/index
      std::stringstream message;
      message << "ChemistryPK::SetupAuxiliaryOutput() : "
              << "Output was requested for '" << aux_names_.at(i) 
              << "' (" << name 
              << ") but no chemistry varibles of this name were found.\n";
      Exceptions::amanzi_throw(ChemistryInvalidInput(message.str()));
    } else {
      aux_index_.push_back(index);
    }
  }

  // create the Epetra_MultiVector that will hold the data
  if (nvars > 0) {
    aux_data_ =
        Teuchos::rcp(new Epetra_MultiVector(
            chemistry_state_->mesh_maps()->cell_map(false), nvars));
  } else {
    aux_data_ = Teuchos::null;
  }
}  // end SetupAuxiliaryOutput()

void Alquimia_Chemistry_PK::CopyAmanziStateToAlquimia(
    const int cell_id,
    Teuchos::RCP<const Epetra_MultiVector> aqueous_components) {

  AlquimiaState* alquimia_state = &chem_data_.state;

  alquimia_state->water_density = (*chemistry_state_->water_density())[cell_id];
  alquimia_state->porosity = (*chemistry_state_->porosity())[cell_id];

  // Alquimia breaks the total aqueous component concentrations into mobile 
  // and immobile parts. The "total sorbed" components are immobile, and 
  // the remainder of the total are mobile.
  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_components = (*aqueous_components)[c];
    double total = cell_components[cell_id];
    double immobile = 0.0;
    if (using_sorption()) {
      double* cell_total_sorbed = (*chemistry_state_->total_sorbed())[c];
      immobile = cell_total_sorbed[cell_id];
    }  // end if(using_sorption)
    double mobile = total - immobile;
    alquimia_state->total_mobile.data[c] = mobile;
    alquimia_state->total_immobile.data[c] = immobile;
  }

  //
  // minerals
  //
  assert(alquimia_state.mineral_volume_fraction.size == number_minerals());
  assert(alquimia_state.mineral_specific_surface_area.size == number_minerals());
  for (unsigned int m = 0; m < number_minerals(); m++) {
    double* cell_minerals = (*chemistry_state_->mineral_volume_fractions())[m];
    alquimia_state->mineral_volume_fraction.data[m] = cell_minerals[cell_id];
    if (chemistry_state_->mineral_specific_surface_area() != Teuchos::null) {
      double* cells_ssa = (*chemistry_state_->mineral_specific_surface_area())[m];
      alquimia_state->mineral_specific_surface_area.data[m] = cells_ssa[cell_id];
    }
  }

  //
  // ion exchange
  //
  assert(alquimia_state->cation_exchange_capacity.size == number_ion_exchange_sites());
  if (number_ion_exchange_sites() > 0) {
    for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) {
      double* cell_ion_exchange_sites = (*chemistry_state_->ion_exchange_sites())[i];
      alquimia_state->cation_exchange_capacity[i] = cell_ion_exchange_sites[cell_id];
    }
  }
  
  //
  // surface complexation
  //
  if (number_sorption_sites() > 0) {
    assert(number_sorption_sites() == alquimia_state->surface_site_density.size);
    for (int s = 0; s < number_sorption_sites(); ++s) {
      // FIXME: Need site density names, too?
      double* cell_sorption_sites = (*chemistry_state_->sorption_sites())[s];
      alquimia_state->surface_site_density.data[s] = cell_sorption_sites[cell_id];
      // TODO(bandre): need to save surface complexation free site conc here!
    }
  }

  // Free ion concentrations, activity coefficients, ion exchange 
  // ref cation sorbed concentrations, and surface complexation free 
  // site concentrations are stored as auxiliary variables internally within
  // Alquimia, but are stored by Amanzi and used as initial estimates for the 
  // new equilibrium conditions as needed. For now, we assume that the data 
  // is ordered as it is in the PFlotran backend, which is as follows:
  //
  // doubles:
  //   free ion conc <N_primary>
  //   primary_activity_coeff <N_primary>
  //   secondary_activity_coeff <N_aqueous_complexes>
  //   ion exchange ref cation sorbed conc <N_ion_exchange_sites>
  //   surface complexation free site conc <N_surface_sites>

  unsigned int offset = 0;
  //
  // free ion concentrations
  //
  for (unsigned int c = 0; c < number_aqueous_components(); c++) { // FIXME: Is this the right number?
    double* cell_free_ion = (*chemistry_state_->free_ion_species())[c];
    chem_data_->aux_data.aux_doubles.data[offset] = cell_free-ion[cell_id];
  }

  //
  // activity coefficients
  //
  for (int i = 0; i < chemistry_state_->primary_activity_coeff()->length(); ++i, ++offset) {
    double* cells = (*chemistry_state_->primary_activity_coeff())[i];
    chem_data_->aux_data.aux_doubles.data[offset] = cells[cell_id];
  }
  for (int i = 0; i < chemistry_state_->secondary_activity_coeff()->length(); ++i, ++offset) {
    double* cells = (*chemistry_state_->secondary_activity_coeff())[i];
    chem_data_->aux_data.aux_doubles.data[offset] = cells[cell_id];
  }

  //
  // ion exchange ref cation concentrations
  //
  for (int i = 0; i < chemistry_state_->ion_exchange_ref_cation_conc()->length(); ++i, ++offset) {
    double* cells = (*chemistry_state_->ion_exchange_ref_cation_conc())[i];
    chem_data_->aux_data.aux_doubles.data[offset] = cells[cell_id];
  }

  //
  // surface complexation
  // 
  for (int i = 0; i < chemistry_state_->surface_complex_free_site_conc()->length(); ++i, ++offset) {
    double* cells = (*chemistry_state_->surface_complex_free_site_conc())[i];
    chem_data_->aux_data.aux_doubles.data[offset] = cells[cell_id];
  }

}  // end CopyAmanziStateToAlquimia()

void Alquimia_Chemistry_PK::CopyAmanziMaterialPropertiesToAlquimiaMaterialProperties(
    const int cell_id,
    Teuchos::RCP<const Epetra_MultiVector> aqueous_components) {

  AlquimiaMaterialProperties* alquimia_mat_props = &chem_data_.material_properties;

  alquimia_mat_props->volume = (*chemistry_state_->volume())[cell_id];
  alquimia_mat_props->saturation = (*chemistry_state_->water_saturation())[cell_id];

  //
  // sorption isotherms
  //
  if (using_sorption_isotherms()) {
    for (unsigned int i = 0; i < number_aqueous_components(); ++i) {
      double* cell_data = (*chemistry_state_->isotherm_kd())[i];
      alquimia_mat_props->isotherm_kd.data[i] = cell_data[cell_id];
      
      cell_data = (*chemistry_state_->isotherm_freundlich_n())[i];
      alquimia_mat_props->isotherm_freundlich_n[i] = cell_data[cell_id];
      
      cell_data = (*chemistry_state_->isotherm_langmuir_b())[i];
      alquimia_mat_props->isotherm_langmuir_b[i] = cell_data[cell_id];
    }
  }

}  // end CopyAmanziMaterialPropertiesToAlquimiaMaterialProperties()

void CopyAmanziGeochemicalConditionsToAlquimia(const int cell_id,
                                               Teuchos::RCP<const Epetra_MultiVector> aqueous_components)
{
  // FIXME: Need to understand how Amanzi conveys geochemical conditions before 
  // FIXME: proceeding with this.
}

void Alquimia_Chemistry_PK::CopyAlquimiaStateToAmanzi(const int cell_id) {
  AlquimiaState* alquimia_state = &chem_data_.state;

  // If the chemistry has modified the porosity and/or density, it needs to 
  // be updated here.
  (*chemistry_state_->water_density())[cell_id] = alquimia_state->water_density;
  (*chemistry_state_->porosity())[cell_id] = alquimia_state->porosity;

  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double mobile = alquimia_state->total_mobile.data[c];
    double immobile = alquimia_state->total_mobile.data[c];

    // NOTE: We place our mobile concentrations into the chemistry state's
    // NOTE: total component concentrations, not the aqueous concentrations.
    // NOTE: I assume this has to do with our operator splitting technique.
    double total = mobile + immobile;
    double* cell_components = (*chemistry_state_->total_component_concentration())[c];
    cell_components[cell_id] = total;

    if (using_sorption())
    {
      double* cell_total_sorbed = (*chemistry_state_->total_sorbed())[c];
      cell_total_sorbed[cell_id] = immobile;
    }
  }

  //
  // minerals
  //
  for (unsigned int m = 0; m < number_minerals(); m++) {
    double* cell_minerals = (*chemistry_state_->mineral_volume_fractions())[m];
    cell_minerals[cell_id] = alquimia_state->mineral_volume_fraction.data[m];
    if (chemistry_state_->mineral_specific_surface_area() != Teuchos::null) {
      cell_minerals = (*chemistry_state_->mineral_specific_surface_area())[m];
      cell_minerals[cell_id] = alquimia_state->mineral_specific_surface_area.data[m];
    }
  }

  //
  // ion exchange
  //
  for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) {
    double* cell_ion_exchange_sites = (*chemistry_state_->ion_exchange_sites())[i];
    cell_ion_exchange_sites[cell_id] = alquimia_state->cation_exchange_capacity.data[i];
  }

  //
  // surface complexation
  //
  for (unsigned int i = 0; i < number_sorption_sites(); i++) {
    double* cell_sorption_sites = (*chemistry_state_->sorption_sites())[i];
    cell_sorption_sites[cell_id] = alquimia_state->surface_site_density.data[i];
  }

  // Free ion concentrations, activity coefficients, ion exchange 
  // ref cation sorbed concentrations, and surface complexation free 
  // site concentrations are all computed internally by Alquimia.
  // We need to extract them from Alquimia's auxiliary arrays.
  // For now, we assume that the data is ordered as it is in the 
  // PFlotran backend, which is as follows:
  //
  // doubles:
  //   free ion conc <N_primary>
  //   primary_activity_coeff <N_primary>
  //   secondary_activity_coeff <N_aqueous_complexes>
  //   ion exchange ref cation sorbed conc <N_ion_exchange_sites>
  //   surface complexation free site conc <N_surface_sites>

  //
  // free ion concentrations
  //
  for (unsigned int c = 0; c < number_aqueous_components(); c++) { // FIXME: Is this the right number?
    double* cell_free_ion = (*chemistry_state_->free_ion_species())[c];
    cell_free_ion[cell_id] = chem_data_->aux_data.aux_doubles.data[offset];
  }

  //
  // activity coefficients
  //
  for (int i = 0; i < chemistry_state_->primary_activity_coeff()->length(); ++i, ++offset) {
    double* cells = (*chemistry_state_->primary_activity_coeff())[i];
    cells[cell_id] = chem_data_->aux_data.aux_doubles.data[offset];
  }
  for (int i = 0; i < chemistry_state_->secondary_activity_coeff()->length(); ++i, ++offset) {
    double* cells = (*chemistry_state_->secondary_activity_coeff())[i];
    cells[cell_id] = chem_data_->aux_data.aux_doubles.data[offset];
  }

  //
  // ion exchange ref cation concentrations
  //
  for (int i = 0; i < chemistry_state_->ion_exchange_ref_cation_conc()->length(); ++i, ++offset) {
    double* cells = (*chemistry_state_->ion_exchange_ref_cation_conc())[i];
    cells[cell_id] = chem_data_->aux_data.aux_doubles.data[offset];
  }

  //
  // surface complexation
  // 
  for (int i = 0; i < chemistry_state_->surface_complex_free_site_conc()->length(); ++i, ++offset) {
    double* cells = (*chemistry_state_->surface_complex_free_site_conc())[i];
    cells[cell_id] = chem_data_->aux_data.aux_doubles.data[offset];
  }

  //
  // Auxiliary output.
  //
  if (aux_data_ != Teuchos::nul) {
    for (unsigned int i = 0; i < aux_names_.size(); i++) {
      if (aux_names_.at(i) == "pH") {
        double* cell_aux_data = (*aux_data_)[i];
        cell_aux_data[cell_id] = chem_data_.aux_data.pH;
      }
      // FIXME: For now, we don't support mineral saturation index or 
      // FIXME: mineral reaction rate, but when we do, these need to 
      // FIXME: be updated here as well.
    }
  }

}  // end CopyAlquimiaStateToAmanzi()

void Alquimia_Chemistry_PK::CopyAlquimiaMaterialPropertiesToAmanzi(const int cell_id) {

  AlquimiaMaterialProperties* alquimia_mat_props = &chem_data_.material_properties;

  // FIXME: Do we need to copy over volume and saturation, too?

  //
  // sorption isotherms
  //
  if (using_sorption_isotherms()) {
    for (unsigned int i = 0; i < number_aqueous_components(); ++i) {
      double* cell_data = (*chemistry_state_->isotherm_kd())[i];
      cell_data[cell_id] = alquimia_mat_props->isotherm_kd.data[i];

      cell_data = (*chemistry_state_->isotherm_freundlich_n())[i];
      cell_data[cell_id] = alquimia_mat_props->isotherm_freundlich_n.data[i];

      cell_data = (*chemistry_state_->isotherm_langmuir_b())[i];
      cell_data[cell_id] = alquimia_mat_props->isotherm_langmuir_b.data[i];
    }
  }
}

/*******************************************************************************
 **
 **  MPC interface functions
 **
 ******************************************************************************/
Teuchos::RCP<Epetra_MultiVector> Alquimia_Chemistry_PK::get_total_component_concentration(void) const {
  return chemistry_state_->total_component_concentration();
}  // end get_total_component_concentration()


// This helper advances the solution on a single cell within Amanzi's state.
// It returns the number of iterations taken to obtain the advanced solution, 
// or -1 if an error occurred.
int Alquimia_Chemistry_PK::AdvanceSingleCell(int cellIndex) {
  // create a backup copy of the components
  //      chem_->CopyComponents(beaker_components_, &beaker_components_copy_);
  // chemistry computations for this cell
  chem.ReactionStepOperatorSplit(&chem_data_.engine_state,
                                 &delta_time,
                                 &chem_data_.material_properties,
                                 &chem_data_.state,
                                 &chem_data_.aux_data,
                                 &chem_status_);
  if (!chem_status_converged)
  {
    return -1;
  }
  else
  {
    int num_iterations = chem_status_.num_newton_iterations;

    // Retrieve the auxiliary output data.
    chem_.GetAuxiliaryOutput(&chem_data_.engine_state,
        &chem_data_.material_properties,
        &chem_data_.state,
        &chem_data_.aux_data,
        &chem_data_.aux_output,
        &chem_status_);

    return num_iterations;
  }
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
    Teuchos::RCP<const Epetra_MultiVector> total_component_concentration_star) {
  if (debug()) {
    chem_out->Write(kVerbose, "  Alquimia_Chemistry_PK::advance() : advancing the chemistry process model...\n");
  }

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

  int ierr(0);
  for (int cell = 0; cell < num_cells; cell++) {
    int num_iterations = AdvanceSingleCell(cell);
    if (num_iterations > 0)
    {
      //std::stringstream message;
      //message << "--- " << cell << "\n";
      //beaker_components_.Display(message.str().c_str());
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
    else
    {
      ierr = 1;
    }

    // update this cell's data in the arrays
    if (ierr == 0) CopyAlquimiaStateToAmanzi(cell);

    // TODO(bandre): was porosity etc changed? copy someplace

#ifdef GLENN_DEBUG
    if (cell % (num_cells / 10) == 0) {
      std::cout << "  " << cell * 100 / num_cells
                << "%" << std::endl;
    }
#endif
  }  // for(cells)
  int recv(0);
  chemistry_state_->mesh_maps()->get_comm()->MaxAll(&ierr, &recv, 1);
  if (recv != 0) {
    AmanziException geochem_error("Error in Alquimia_Chemistry_PK::advance");
    Exceptions::amanzi_throw(geochem_error); 
  }  
  
#ifdef GLENN_DEBUG
  if (debug() == kDebugChemistryProcessKernel) {
    std::stringstream message;
    message << "  Alquimia_Chemistry_PK::advance() : "
            << "max iterations - " << max_iterations << " " << "  cell id: " << imax << std::endl;
    message << "  Alquimia_Chemistry_PK::advance() : "
            << "min iterations - " << min_iterations << " " << "  cell id: " << imin << std::endl;
    message << "  Alquimia_Chemistry_PK::advance() : "
            << "ave iterations - " << static_cast<float>(ave_iterations) / num_cells << std::endl;
  }
#endif

#if 0
  if (debug() == kDebugChemistryProcessKernel) {
    // dumping the values of the final cell. not very helpful by itself,
    // but can be move up into the loops....
    chem_->DisplayTotalColumnHeaders(display_free_columns_);
    chem_->DisplayTotalColumns(current_time_, beaker_components_, true);
  }
#endif
}  // end advance()


// the MPC will call this function to signal to the
// process kernel that it has accepted the
// state update, thus, the PK should update
// possible auxilary state variables here
void Alquimia_Chemistry_PK::commit_state(Teuchos::RCP<Chemistry_State> chem_state,
                                const double& delta_time) {
  if (debug() == kDebugChemistryProcessKernel) {
    chem_out->Write(kVerbose,
                    "  Alquimia_Chemistry_PK::commit_state() : Committing internal state.\n");
  }

  saved_time_ += delta_time;

  // FIXME: I'm not sure there's anything else we have to do here.
#if 0
  if (debug() && false) {
    chem_->Speciate(&beaker_components_, beaker_parameters_);
    chem_->DisplayResults();
    chem_->DisplayTotalColumnHeaders(display_free_columns_);
    chem_->DisplayTotalColumns(saved_time_, beaker_components_, true);
  }
#endif
}  // end commit_state()



Teuchos::RCP<Epetra_MultiVector> Alquimia_Chemistry_PK::get_extra_chemistry_output_data() {
  // This vector is updated during the initialization and advance of 
  // the geochemistry, so we simply return it here.
  return aux_data_;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
