/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "chemistry_pk.hh"

#include <string>
#include <algorithm>

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

#include "chemistry_state.hh"
#include "simple_thermo_database.hh"
#include "beaker.hh"
#include "chemistry_output.hh"
#include "chemistry_verbosity.hh"
#include "chemistry_exception.hh"

#include "Mesh.hh"
#include "State.hpp"
#include "errors.hh"
#include "exceptions.hh"

namespace amanzi {
namespace chemistry {

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

// global ChemistryOutput object in the amanzi::chemisry
// namespace that will be used by the chemistry library
extern ChemistryOutput* chem_out;

Chemistry_PK::Chemistry_PK(const Teuchos::ParameterList& param_list,
                           Teuchos::RCP<Chemistry_State> chem_state)
    : debug_(false),
      display_free_columns_(false),
      max_time_step_(9.9e9),
      chemistry_state_(chem_state),
      parameter_list_(param_list),
      chem_(NULL),
      current_time_(0.0),
      saved_time_(0.0) {

  // NOTE: we want the chemistry_pk and chem lib verbose by default so
  // we can debug them. add the "silent" level to the input file or
  // have the MPC do an if (my_mpi_process) to prevent generating a
  // lot of output in parallel machines.
  SetupDefaultChemistryOutput();
  chem_out->AddLevel("verbose");
}  // end Chemistry_PK()

Chemistry_PK::~Chemistry_PK() {
  delete chem_;
  delete chem_out;
}  // end ~Chemistry_PK()

void Chemistry_PK::InitializeChemistry(void) {
  if (debug()) {
    std::cout << "  Chemistry_PK::InitializeChemistry()" << std::endl;
  }

  XMLParameters();

  // modify output levels based on xml input
  // for (levels) {chem_out->AddLevel()}
  // debugging
  if (false) {
    chem_out->RemoveLevel("silent");
    chem_out->AddLevel("some debug level");
    chem_out->set_use_stdout(true);
  }

  // TODO/NOTE(bandre): Turning off warnings to silence charge balance
  // warnings. Need a better way to do this...
  chem_out->RemoveLevel("warning");

  // TODO: some sort of check of the state object to see if mineral_ssa,
  // CEC, site density, etc is present.

  // initial conditions for minerals etc should be handled by the
  // state/chemistry_state object before we reach this point. We just
  // resize our local memory for migrating data here.

  SizeBeakerStructures();

  // copy the first cell data into the beaker storage for
  // initialization purposes
  int cell = 0;
  CopyCellStateToBeakerStructures(
      cell, 
      chemistry_state_->total_component_concentration());

  // finish setting up & testing the chemistry object
  try {
    chem_->set_debug(false);
    chem_out->Write(kVerbose, "ChemistryPK: Initializing chemistry in cell 0...\n");
    chem_->Setup(beaker_components_, beaker_parameters_);
    chem_->Display();
    // solve for initial free-ion concentrations
    chem_out->Write(kVerbose, "ChemistryPK: Initial speciation calculations in cell 0...\n");
    chem_->Speciate(&beaker_components_, beaker_parameters_);
    if (debug()) {
      chem_out->Write(kVerbose, "\nTest solution of initial conditions in cell 0:\n");
      chem_->DisplayResults();
    }
  } catch (ChemistryException& geochem_error) {
    chem_->DisplayResults();
    std::cout << geochem_error.what() << std::endl;
    Exceptions::amanzi_throw(geochem_error);
  }

  // TODO(bandre): at this point we should know about any additional
  // storage that chemistry needs...
  chemistry_state_->AllocateAdditionalChemistryStorage(beaker_components_);

  SetupAuxiliaryOutput();

  // now loop through all the cells and initialize
  chem_out->Write(kVerbose, "ChemistryPK: Initializing chemistry in all cells...\n");
  int num_cells = chemistry_state_->porosity()->MyLength();
  for (cell = 0; cell < num_cells; ++cell) {
    if (debug()) {
      std::cout << "Initial speciation in cell " << cell << std::endl;
    }

    CopyCellStateToBeakerStructures(
        cell,
        chemistry_state_->total_component_concentration());

    try {
      //chem_->DisplayTotalColumns(static_cast<double>(-cell), beaker_components_, display_free_columns_);
      // solve for initial free-ion concentrations
      chem_->Speciate(&beaker_components_, beaker_parameters_);
      CopyBeakerStructuresToCellState(cell);
    } catch (ChemistryException& geochem_error) {
      std::cout << "ChemistryPK::InitializeChemistry(): " 
                << "An error occured while initializing chemistry in cell: " 
                << cell << ".\n" << geochem_error.what() << std::endl;
      Exceptions::amanzi_throw(geochem_error);
    }

#ifdef GLENN_DEBUG
    if (cell % (num_cells / 10) == 0) {
      std::cout << "  " << cell * 100 / num_cells
                << "%" << std::endl;
    }
#endif
  }
  chem_out->Write(kVerbose, "ChemistryPK::InitializeChemistry(): initialization was successful.\n");
}  // end InitializeChemistry()

/*******************************************************************************
 **
 **  initialization helper functions
 **
 ******************************************************************************/
void Chemistry_PK::XMLParameters(void) {
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
  // determine the format of the database file then create the database!
  //
  //--------------------------------------------------------------------------
  if (parameter_list_.isParameter("Thermodynamic Database Format")) {
    std::string database_format =
        parameter_list_.get<std::string>("Thermodynamic Database Format");
    // create the appropriate chemistry object
    if (database_format == "simple") {
      chem_ = new SimpleThermoDatabase();
    } else {
      // invalid database format...
      std::ostringstream error_stream;
      error_stream << ChemistryException::kChemistryError;
      error_stream << "Chemistry_PK::XMLParameters(): \n";
      error_stream << "  Input parameter 'Thermodynamic Database Format' must be 'simple'.\n";
      Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
    }
  } else {
    std::ostringstream error_stream;
    error_stream << ChemistryException::kChemistryError;
    error_stream << "Chemistry_PK::XMLParameters(): \n";
    error_stream << "  Input parameter 'Thermodynamic Database Format' must be specified.\n";
    Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
  }

  beaker_parameters_ = chem_->GetDefaultParameters();

  //--------------------------------------------------------------------------
  //
  // thermo file name
  //
  //--------------------------------------------------------------------------
  if (parameter_list_.isParameter("Thermodynamic Database File")) {
    beaker_parameters_.thermo_database_file =
        parameter_list_.get<std::string>("Thermodynamic Database File");
  } else {
    std::ostringstream error_stream;
    error_stream << ChemistryException::kChemistryError;
    error_stream << "Chemistry_PK::XMLParameters(): \n";
    error_stream << "  Input parameter 'Thermodynamic Database File' must be specified.\n";
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
      error_stream << ChemistryException::kChemistryError;
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


void Chemistry_PK::SetupAuxiliaryOutput(void) {
  // requires that Beaker::Setup() has already been called!
  if (debug()) {
    std::cout << "  Chemistry_PK::SetupAuxiliaryOutput()" << std::endl;
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


void Chemistry_PK::SizeBeakerStructures(void) {
  // initialize the beaker component data structure

  // NOTE: The beaker already has data for site density, sorption
  // isotherms, ssa. If we want to use that single global value, then
  // we leave these arrays empty as a flag to the beaker to use its
  // own value. If we want to over ride the global chemistry value
  // with cell by cell data, then we resize the containers here.

  beaker_components_.total.resize(number_aqueous_components(), 0.0);
  beaker_components_.free_ion.resize(number_aqueous_components(), 1.0e-9);
  beaker_components_.mineral_volume_fraction.resize(number_minerals(), 0.0);

  if (using_sorption()) {
    beaker_components_.total_sorbed.resize(number_total_sorbed(), 0.0);
  } else {
    beaker_components_.total_sorbed.clear();
  }

  if (number_minerals() > 0) {
    beaker_components_.mineral_specific_surface_area.resize(number_minerals(), 0.0);
  } else {
    beaker_components_.mineral_specific_surface_area.clear();
  }

  if (number_ion_exchange_sites() > 0) {
    beaker_components_.ion_exchange_sites.resize(number_ion_exchange_sites(), 0.0);
  } else {
    beaker_components_.ion_exchange_sites.clear();
  }

  if (number_sorption_sites() > 0) {
    beaker_components_.surface_site_density.resize(number_sorption_sites(), 0.0);
  } else {
    beaker_components_.surface_site_density.clear();
  }

  if (using_sorption_isotherms()) {
    beaker_components_.isotherm_kd.resize(number_aqueous_components(), 0.0);
    beaker_components_.isotherm_freundlich_n.resize(number_aqueous_components(), 0.0);
    beaker_components_.isotherm_langmuir_b.resize(number_aqueous_components(), 0.0);
  } else {
    beaker_components_.isotherm_kd.clear();
    beaker_components_.isotherm_freundlich_n.clear();
    beaker_components_.isotherm_langmuir_b.clear();
  }
}  // end SizeBeakerStructures()


void Chemistry_PK::CopyCellStateToBeakerStructures(
    const int cell_id,
    Teuchos::RCP<const Epetra_MultiVector> aqueous_components) {
  // NOTE: want the aqueous totals value calculated from transport
  // (aqueous_components), not the value stored in state!

  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_components = (*aqueous_components)[c];
    beaker_components_.total.at(c) = cell_components[cell_id];
  }

  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_free_ion = (*chemistry_state_->free_ion_species())[c];
    beaker_components_.free_ion.at(c) = cell_free_ion[cell_id];
  }

  //
  // activity coefficients
  //
  for (unsigned int i = 0; i < beaker_components_.primary_activity_coeff.size(); ++i) {
    double* cells = (*chemistry_state_->primary_activity_coeff())[i];
    beaker_components_.primary_activity_coeff.at(i) = cells[cell_id];
  }

  for (unsigned int i = 0; i < beaker_components_.secondary_activity_coeff.size(); ++i) {
    double* cells = (*chemistry_state_->secondary_activity_coeff())[i];
    beaker_components_.secondary_activity_coeff.at(i) = cells[cell_id];
  }

  //
  // minerals
  //
  for (unsigned int m = 0; m < number_minerals(); m++) {
    double* cell_minerals = (*chemistry_state_->mineral_volume_fractions())[m];
    beaker_components_.mineral_volume_fraction[m] = cell_minerals[cell_id];
    if (chemistry_state_->mineral_specific_surface_area() != Teuchos::null) {
      double* cells_ssa = (*chemistry_state_->mineral_specific_surface_area())[m];
      beaker_components_.mineral_specific_surface_area.at(m) = cells_ssa[cell_id];
    }
  }

  //
  // general sorption storage
  //
  if (using_sorption()) {
    for (unsigned int c = 0; c < number_aqueous_components(); c++) {
      double* cell_total_sorbed = (*chemistry_state_->total_sorbed())[c];
      beaker_components_.total_sorbed.at(c) = cell_total_sorbed[cell_id];
    }
  }  // end if(using_sorption)

  //
  // ion exchange
  //
  if (number_ion_exchange_sites() > 0) {
    // TODO: only allow one ion exchange site at the moment!
    for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) {
      double* cell_ion_exchange_sites = (*chemistry_state_->ion_exchange_sites())[i];
      beaker_components_.ion_exchange_sites[i] = cell_ion_exchange_sites[cell_id];
      // TODO(bandre): need to save ion exchange ref cation conc here!
    }
  }
  
  for (unsigned int i = 0; i < beaker_components_.ion_exchange_ref_cation_conc.size(); ++i) {
    double* cells = (*chemistry_state_->ion_exchange_ref_cation_conc())[i];
    beaker_components_.ion_exchange_ref_cation_conc.at(i) = cells[cell_id];
  }

  //
  // surface complexation
  //
  if (number_sorption_sites() > 0) {
    for (int s = 0; s < number_sorption_sites(); ++s) {
      double* cell_sorption_sites = 
          (*chemistry_state_->sorption_sites())[s];
      beaker_components_.surface_site_density[s] = cell_sorption_sites[cell_id];
      // TODO(bandre): need to save surface complexation free site conc here!
    }
  }

  for (unsigned int i = 0; i < beaker_components_.surface_complex_free_site_conc.size(); ++i) {
    double* cells = (*chemistry_state_->surface_complex_free_site_conc())[i];
    beaker_components_.surface_complex_free_site_conc.at(i) = cells[cell_id];
  }

  //
  // sorption isotherms
  //
  if (using_sorption_isotherms()) {
    for (unsigned int i = 0; i < number_aqueous_components(); ++i) {
      double* cell_data = (*chemistry_state_->isotherm_kd())[i];
      beaker_components_.isotherm_kd.at(i) = cell_data[cell_id];
      
      cell_data = (*chemistry_state_->isotherm_freundlich_n())[i];
      beaker_components_.isotherm_freundlich_n.at(i) = cell_data[cell_id];
      
      cell_data = (*chemistry_state_->isotherm_langmuir_b())[i];
      beaker_components_.isotherm_langmuir_b.at(i) = cell_data[cell_id];
    }
  }

  // copy data from state arrays into the beaker parameters
  beaker_parameters_.water_density = (*chemistry_state_->water_density())[cell_id];
  beaker_parameters_.porosity = (*chemistry_state_->porosity())[cell_id];
  beaker_parameters_.saturation = (*chemistry_state_->water_saturation())[cell_id];
  beaker_parameters_.volume = (*chemistry_state_->volume())[cell_id];

}  // end CopyCellStateToBeakerStructures()


void Chemistry_PK::CopyBeakerStructuresToCellState(const int cell_id) {
  // copy data from the beaker back into the state arrays

  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_components = (*chemistry_state_->total_component_concentration())[c];
    cell_components[cell_id] = beaker_components_.total.at(c);
  }

  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_free_ion = (*chemistry_state_->free_ion_species())[c];
    cell_free_ion[cell_id] = beaker_components_.free_ion.at(c);
  }

  //
  // activity coefficients
  //
  for (unsigned int i = 0; i < beaker_components_.primary_activity_coeff.size(); ++i) {
      double* cells = (*chemistry_state_->primary_activity_coeff())[i];
      cells[cell_id] = beaker_components_.primary_activity_coeff.at(i);
  }
  for (unsigned int i = 0; i < beaker_components_.secondary_activity_coeff.size(); ++i) {
    double* cells = (*chemistry_state_->secondary_activity_coeff())[i];
    cells[cell_id] =  beaker_components_.secondary_activity_coeff.at(i);
  }

  //
  // minerals
  //
  for (unsigned int m = 0; m < number_minerals(); m++) {
    double* cell_minerals = (*chemistry_state_->mineral_volume_fractions())[m];
    cell_minerals[cell_id] = beaker_components_.mineral_volume_fraction.at(m);
    if (chemistry_state_->mineral_specific_surface_area() != Teuchos::null) {
      cell_minerals = (*chemistry_state_->mineral_specific_surface_area())[m];
      cell_minerals[cell_id] = beaker_components_.mineral_specific_surface_area.at(m);
    }
  }

  //
  // sorption
  //
  if (using_sorption()) {
    for (unsigned int c = 0; c < number_aqueous_components(); c++) {
      double* cell_total_sorbed = (*chemistry_state_->total_sorbed())[c];
      cell_total_sorbed[cell_id] = beaker_components_.total_sorbed.at(c);
    }
  }

  //
  // surface complexation
  //
  for (unsigned int i = 0; i < number_sorption_sites(); i++) {
    double* cell_sorption_sites = (*chemistry_state_->sorption_sites())[i];
    cell_sorption_sites[cell_id] = beaker_components_.surface_site_density.at(i);
    // TODO(bandre): need to save surface complexation free site conc here!
  }

  for (unsigned int i = 0; i < beaker_components_.surface_complex_free_site_conc.size(); ++i) {
    double* cells = (*chemistry_state_->surface_complex_free_site_conc())[i];
    cells[cell_id] = beaker_components_.surface_complex_free_site_conc.at(i);
  }

  //
  // ion exchange
  //
  for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) {
    double* cell_ion_exchange_sites = (*chemistry_state_->ion_exchange_sites())[i];
    cell_ion_exchange_sites[cell_id] = beaker_components_.ion_exchange_sites.at(i);
    // TODO(bandre): need to save ion exchange ref cation conc here!
  }

  for (unsigned int i = 0; i < beaker_components_.ion_exchange_ref_cation_conc.size(); ++i) {
    double* cells = (*chemistry_state_->ion_exchange_ref_cation_conc())[i];
    cells[cell_id] = beaker_components_.ion_exchange_ref_cation_conc.at(i);
  }

  //
  // sorption isotherms
  //
  if (using_sorption_isotherms()) {
    for (unsigned int i = 0; i < number_aqueous_components(); ++i) {
      double* cell_data = (*chemistry_state_->isotherm_kd())[i];
      cell_data[cell_id] = beaker_components_.isotherm_kd.at(i);

      cell_data = (*chemistry_state_->isotherm_freundlich_n())[i];
      cell_data[cell_id] = beaker_components_.isotherm_freundlich_n.at(i);

      cell_data = (*chemistry_state_->isotherm_langmuir_b())[i];
      cell_data[cell_id] = beaker_components_.isotherm_langmuir_b.at(i);
    }
  }

  // TODO(bandre): if chemistry can modify the porosity or density,
  // then they should be updated here!

}  // end CopyBeakerStructureToCellState()





/*******************************************************************************
 **
 **  MPC interface functions
 **
 ******************************************************************************/
Teuchos::RCP<Epetra_MultiVector> Chemistry_PK::get_total_component_concentration(void) const {
  return chemistry_state_->total_component_concentration();
}  // end get_total_component_concentration()


/*******************************************************************************
 **
 ** Chemistry_PK::advance()
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

void Chemistry_PK::advance(
    const double& delta_time,
    Teuchos::RCP<const Epetra_MultiVector> total_component_concentration_star) {
  if (debug()) {
    chem_out->Write(kVerbose, "  Chemistry_PK::advance() : advancing the chemistry process model...\n");
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

  for (int cell = 0; cell < num_cells; cell++) {
    // copy state data into the beaker data structures
    CopyCellStateToBeakerStructures(cell, tcc_star);
    try {
      // create a backup copy of the components
      chem_->CopyComponents(beaker_components_, &beaker_components_copy_);
      // chemistry computations for this cell
      int num_iterations = chem_->ReactionStep(&beaker_components_,
                                               beaker_parameters_, delta_time);
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
    } catch (ChemistryException& geochem_error) {
      std::cout << "ERROR: Chemistry_PR::advance() "
                << "cell[" << cell << "]: " << std::endl;
      std::cout << geochem_error.what();
      chem_->DisplayTotalColumnHeaders(display_free_columns_);
      std::cout << "\nFailed Solution" << std::endl;
      std::cout << "  Total Component Concentrations" << std::endl;
      chem_->DisplayTotalColumns(current_time_, beaker_components_, display_free_columns_);
      std::cout << "\nPrevious Solution" << std::endl;
      std::cout << "  Total Component Concentrations" << std::endl;
      chem_->DisplayTotalColumns(current_time_, beaker_components_copy_, display_free_columns_);
      std::cout << std::endl;
      Exceptions::amanzi_throw(geochem_error);
    }

    // update this cell's data in the arrays
    CopyBeakerStructuresToCellState(cell);

    // TODO(bandre): was porosity etc changed? copy someplace

#ifdef GLENN_DEBUG
    if (cell % (num_cells / 10) == 0) {
      std::cout << "  " << cell * 100 / num_cells
                << "%" << std::endl;
    }
#endif
  }  // for(cells)

#ifdef GLENN_DEBUG
  if (debug() == kDebugChemistryProcessKernel) {
    std::stringstream message;
    message << "  Chemistry_PK::advance() : "
            << "max iterations - " << max_iterations << " " << "  cell id: " << imax << std::endl;
    message << "  Chemistry_PK::advance() : "
            << "min iterations - " << min_iterations << " " << "  cell id: " << imin << std::endl;
    message << "  Chemistry_PK::advance() : "
            << "ave iterations - " << static_cast<float>(ave_iterations) / num_cells << std::endl;
  }
#endif

  if (debug() == kDebugChemistryProcessKernel) {
    // dumping the values of the final cell. not very helpful by itself,
    // but can be move up into the loops....
    chem_->DisplayTotalColumnHeaders(display_free_columns_);
    chem_->DisplayTotalColumns(current_time_, beaker_components_, true);
  }
}  // end advance()


// the MPC will call this function to signal to the
// process kernel that it has accepted the
// state update, thus, the PK should update
// possible auxilary state variables here
void Chemistry_PK::commit_state(Teuchos::RCP<Chemistry_State> chem_state,
                                const double& delta_time) {
  if (debug() == kDebugChemistryProcessKernel) {
    chem_out->Write(kVerbose,
                    "  Chemistry_PK::commit_state() : Committing internal state.\n");
  }

  saved_time_ += delta_time;

  if (debug() && false) {
    chem_->Speciate(&beaker_components_, beaker_parameters_);
    chem_->DisplayResults();
    chem_->DisplayTotalColumnHeaders(display_free_columns_);
    chem_->DisplayTotalColumns(saved_time_, beaker_components_, true);
  }
}  // end commit_state()



Teuchos::RCP<Epetra_MultiVector> Chemistry_PK::get_extra_chemistry_output_data() {
  if (aux_data_ != Teuchos::null) {
    int num_cells = chemistry_state_->porosity()->MyLength();

    for (int cell = 0; cell < num_cells; cell++) {
      // populate aux_data_ by copying from the appropriate internal storage
      for (unsigned int i = 0; i < aux_names_.size(); i++) {
        if (aux_names_.at(i) == "pH") {
          double* cell_aux_data = (*aux_data_)[i];
          double* cell_free_ion = (*chemistry_state_->free_ion_species())[aux_index_.at(i)];
          double* activity_coeff = (*chemistry_state_->primary_activity_coeff())[aux_index_.at(i)];
          double activity = cell_free_ion[cell] * activity_coeff[cell];
          cell_aux_data[cell] = -std::log10(activity);
        } else {
          // don't support anything else at this time....
        }
      }  // for(aux_names)
    }  // for(cells)

    // return the multi vector
  }
  return aux_data_;
}


void Chemistry_PK::set_chemistry_output_names(std::vector<string>* names) {
  names->clear();
  for (std::vector<string>::const_iterator name = aux_names_.begin();
       name != aux_names_.end(); name++) {
    names->push_back(*name);
  }
}  // end set_chemistry_output_names()

void Chemistry_PK::set_component_names(std::vector<string>* names) {
  // TODO(bandre): this interface should no longer be used because the
  // xml input file should contain this data in a place state can
  // obtain it
  chem_->GetPrimaryNames(names);
}  // end set_component_names()

}  // namespace chemistry
}  // namespace amanzi
