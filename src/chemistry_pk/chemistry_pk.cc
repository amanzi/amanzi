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
#include "verbosity.hh"
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

Chemistry_PK::Chemistry_PK(const Teuchos::ParameterList& param_list,
                           Teuchos::RCP<Chemistry_State> chem_state)
    : verbosity_(kSilent),
      max_time_step_(9.9e9),
      chemistry_state_(chem_state),
      parameter_list_(param_list),
      chem_(NULL),
      current_time_(0.0),
      saved_time_(0.0) {

  // NOTE: we want the chemistry_pk and chem lib silet by default so
  // they don't generate a lot of output in parallel machines. Can
  // modify in the xml file or after verbosity is set in
  // XMLParameters() for debugging.

}  // end Chemistry_PK()

Chemistry_PK::~Chemistry_PK() {
  delete chem_;
}  // end ~Chemistry_PK()

void Chemistry_PK::InitializeChemistry(void) {
  if (verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "  Chemistry_PK::InitializeChemistry()" << std::endl;
  }

  XMLParameters();
  //set_verbosity(kDebugBeaker);

  // TODO: some sort of check of the state object to see if mineral_ssa,
  // CEC, site density, etc is present.
  set_override_database(true);

  // initial conditions for minerals etc should be handled by the
  // state/chemistry_state object before we reach this point. We just
  // resize our local memory for migrating data here.

  SizeBeakerComponents();

  // copy the first cell data into the beaker storage for
  // initialization purposes
  int cell = 0;
  CopyStateToBeakerParameters(cell);
  CopyCellToBeakerComponents(cell, chemistry_state_->total_component_concentration());

  // finish setting up the chemistry object
  try {
    chem_->verbosity(verbosity());
    chem_->Setup(beaker_components_, beaker_parameters_);
    if (verbosity() > kTerse) {
      chem_->Display();
    }

    // solve for initial free-ion concentrations
    chem_->Speciate(beaker_components_, beaker_parameters_);
    chem_->UpdateComponents(&beaker_components_);
    if (verbosity() > kTerse) {
      std::cout << "\nTest solution of initial conditions in cell 0:"
                << std::endl;
      chem_->DisplayResults();
    }
  } catch (ChemistryException& geochem_error) {
    // TODO(bandre): any errer in the constructor is probably fatal.
    // TODO(bandre): write to cout/cerr here or let the catcher handle it?
    chem_->DisplayResults();
    std::cout << geochem_error.what() << std::endl;
    Exceptions::amanzi_throw(geochem_error);
  }

  CopyBeakerComponentsToCell(cell);

  SetupAuxiliaryOutput();

  // now take care of the remainder
  int num_cells = chemistry_state_->porosity()->MyLength();
  for (int icell = 1; icell < num_cells; icell++) {
    CopyStateToBeakerParameters(icell);
    CopyCellToBeakerComponents(icell, chemistry_state_->total_component_concentration());

    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "Reacting in cell " << icell << std::endl;
    }

    try {
      // solve for initial free-ion concentrations
      chem_->Speciate(beaker_components_, beaker_parameters_);
      chem_->UpdateComponents(&beaker_components_);
    } catch (ChemistryException& geochem_error) {
      std::cout << geochem_error.what() << std::endl;
      Exceptions::amanzi_throw(geochem_error);
    }
    // if successful copy back
    CopyBeakerComponentsToCell(icell);

#ifdef GLENN_DEBUG
    if (icell % (num_cells / 10) == 0) {
      std::cout << "  " << icell * 100 / num_cells
                << "%" << std::endl;
    }
#endif
  }
}  // end InitializeChemistry()

/*******************************************************************************
 **
 **  initialization helper functions
 **
 ******************************************************************************/
void Chemistry_PK::XMLParameters(void) {
  // extract parameters from the xml list and set in the parameters
  // structure

  set_verbosity(static_cast<Verbosity>(parameter_list_.get<int>("Verbosity", 0)));

  if (verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "  Chemistry_PK::XMLParameters()" << std::endl;
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
  // misc other chemistry flags
  //
  // --------------------------------------------------------------------------
  set_max_time_step(parameter_list_.get<double>("Max Time Step (s)", 9.9e9));

}  // end XMLParameters()


void Chemistry_PK::SetupAuxiliaryOutput(void) {
  // requires that Beaker::Setup() has already been called!
  if (verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "  Chemistry_PK::SetupAuxiliaryOutput()" << std::endl;
  }
  // TODO(bandre): temporary hard coding of auxillary output names, needs to
  // come from the input file eventually
  aux_names_.clear();
  aux_names_.push_back("pH");

  unsigned int nvars = aux_names_.size();
  std::string name;
  aux_index_.clear();
  for (unsigned int i = 0; i < nvars; i++) {
    if (aux_names_.at(i) == "pH") {
      name = "H+";
    } else {
      name = aux_names_.at(i);
    }
    aux_index_.push_back(chem_->GetPrimaryIndex(name));
    // check to make sure it is not -1, an invalid name/index
  }

  // create the Epetra_MultiVector that will hold the data
  aux_data_ =
      Teuchos::rcp(new Epetra_MultiVector(
          chemistry_state_->mesh_maps()->cell_map(false), nvars));
}  // end SetupAuxiliaryOutput()


void Chemistry_PK::SizeBeakerComponents(void) {
  // initialize the beaker component data structure
  beaker_components_.total.clear();
  beaker_components_.free_ion.clear();
  beaker_components_.minerals.clear();
  beaker_components_.ion_exchange_sites.clear();
  beaker_components_.total_sorbed.clear();

  beaker_components_.total.resize(number_aqueous_components(), 0.0);
  beaker_components_.minerals.resize(number_minerals(), 0.0);
  beaker_components_.free_ion.resize(number_aqueous_components(), 1.0e-9);

  if (using_sorption()) {
    beaker_components_.total_sorbed.resize(number_total_sorbed(), 0.0);
  }

  if (override_database()) {
    beaker_parameters_.mineral_specific_surface_area.resize(number_minerals(), 0.0);
    if (number_ion_exchange_sites() > 0) {
      beaker_components_.ion_exchange_sites.resize(number_ion_exchange_sites(), 0.0);
    }

    if (number_sorption_sites() > 0) {
      beaker_parameters_.sorption_site_density.resize(number_sorption_sites(), 0.0);
    }
  }  // if(override)
}  // end SizeBeakerComponents()


void Chemistry_PK::CopyCellToBeakerComponents(
    const int cell_id,
    Teuchos::RCP<const Epetra_MultiVector> aqueous_components) {
  // copy component data from the cell arrays into beaker component
  // structure

  // Note: total uses the passed in value from transport, not the
  // current storage value
  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_components = (*aqueous_components)[c];
    beaker_components_.total[c] = cell_components[cell_id];
  }

  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_free_ion = (*chemistry_state_->free_ion_species())[c];
    beaker_components_.free_ion[c] = cell_free_ion[cell_id];
  }

  for (unsigned int m = 0; m < number_minerals(); m++) {
    double* cell_minerals = (*chemistry_state_->mineral_volume_fractions())[m];
    beaker_components_.minerals[m] = cell_minerals[cell_id];
  }

  for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) {
    double* cell_ion_exchange_sites = (*chemistry_state_->ion_exchange_sites())[i];
    beaker_components_.ion_exchange_sites[i] = cell_ion_exchange_sites[cell_id];
  }

  if (using_sorption()) {
    for (unsigned int c = 0; c < number_aqueous_components(); c++) {
      double* cell_total_sorbed = (*chemistry_state_->total_sorbed())[c];
      beaker_components_.total_sorbed[c] = cell_total_sorbed[cell_id];
    }
  }  // end if(using_sorption)
}  // end CopyCellToBeakerComponents()


void Chemistry_PK::CopyBeakerComponentsToCell(const int cell_id) {
  // copy data from the beaker back into the state arrays

  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_components = (*chemistry_state_->total_component_concentration())[c];
    cell_components[cell_id] = beaker_components_.total[c];
  }

  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_free_ion = (*chemistry_state_->free_ion_species())[c];
    cell_free_ion[cell_id] = beaker_components_.free_ion[c];
  }

  for (unsigned int m = 0; m < number_minerals(); m++) {
    double* cell_minerals = (*chemistry_state_->mineral_volume_fractions())[m];
    cell_minerals[cell_id] = beaker_components_.minerals[m];
  }

  for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) {
    double* cell_ion_exchange_sites = (*chemistry_state_->ion_exchange_sites())[i];
    cell_ion_exchange_sites[cell_id] = beaker_components_.ion_exchange_sites[i];
  }

  if (using_sorption()) {
    for (unsigned int c = 0; c < number_aqueous_components(); c++) {
      double* cell_total_sorbed = (*chemistry_state_->total_sorbed())[c];
      cell_total_sorbed[cell_id] = beaker_components_.total_sorbed[c];
    }
    if (number_sorption_sites() > 0) {
      for (unsigned int i = 0; i < number_sorption_sites(); i++) {
        double* cell_sorption_sites = (*chemistry_state_->sorption_sites())[i];
        cell_sorption_sites[cell_id] = beaker_components_.total_sorbed[i];
      }
    }
  }  // end if(using_sorption)
}  // end CopyBeakerComponentsToCell()


void Chemistry_PK::CopyStateToBeakerParameters(const int cell_id) {
  // copy data from state arrays into the beaker parameters
  beaker_parameters_.water_density = (*chemistry_state_->water_density())[cell_id];
  beaker_parameters_.porosity = (*chemistry_state_->porosity())[cell_id];
  beaker_parameters_.saturation = (*chemistry_state_->water_saturation())[cell_id];
  beaker_parameters_.volume = (*chemistry_state_->volume())[cell_id];
  beaker_parameters_.override_database = override_database();
  if (override_database()) {
    for (int m = 0; m < number_minerals(); ++m) {
      double* cells_ssa = (*chemistry_state_->mineral_specific_surface_area())[m];
      beaker_parameters_.mineral_specific_surface_area.at(m) = cells_ssa[cell_id];
    }

    for (int s = 0; s < number_sorption_sites(); ++s) {
        double* cell_sorption_sites = 
            (*chemistry_state_->sorption_sites())[s];
        beaker_parameters_.sorption_site_density[s] = cell_sorption_sites[cell_id];
    }
    // TODO: only allow one ion exchange site at the moment!
    // if (number_ion_exchange_sites() > 0) {
    //   double* cell_cec = (*chemistry_state_->ion_exchange_sites())[0];
    //   beaker_parameters_.cation_exchange_capacity = cell_cec[cell_id];
    // }
        
  }
}  // end CopyStateToBeakerParameters()



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
  if (chem_->verbosity() == kDebugChemistryProcessKernel) {
    cout << "  Chemistry_PK::advance() : "
         << "advancing the chemistry process model..." << endl;
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
    CopyCellToBeakerComponents(cell, tcc_star);
    CopyStateToBeakerParameters(cell);

    try {
      // chemistry computations for this cell
      chem_->CopyComponents(beaker_components_, &beaker_components_copy_);
      int num_iterations = chem_->ReactionStep(&beaker_components_,
                                               beaker_parameters_, delta_time);
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
      chem_->DisplayTotalColumnHeaders();
      std::cout << "\nFailed Solution" << std::endl;
      std::cout << "  Total Component Concentrations" << std::endl;
      chem_->DisplayTotalColumns(current_time_, beaker_components_);
      // std::cout << "  Free Ion Concentrations" << std::endl;
      // chem_->DisplayTotalColumns(current_time_, beaker_components_);
      // std::cout << "  Total Sorbed Concentrations" << std::endl;
      // chem_->DisplayTotalColumns(current_time_, beaker_components_);
      // #ifdef GLENN_DEBUG
      // TODO(bandre): these cause an exception if called when the above copy is missing
      std::cout << "\nPrevious Solution" << std::endl;
      std::cout << "  Total Component Concentrations" << std::endl;
      chem_->DisplayTotalColumns(current_time_, beaker_components_copy_);
      // std::cout << "  Free Ion Concentrations" << std::endl;
      // chem_->DisplayTotalColumns(current_time_, beaker_components_copy_);
      // std::cout << "  Total Sorbed Concentrations" << std::endl;
      // chem_->DisplayTotalColumns(current_time_, beaker_components_copy_);
      // #endif
      std::cout << std::endl;
      Exceptions::amanzi_throw(geochem_error);
    }

    // update this cell's data in the arrays
    CopyBeakerComponentsToCell(cell);

    // TODO(bandre): was porosity etc changed? copy someplace

#ifdef GLENN_DEBUG
    if (cell % (num_cells / 10) == 0) {
      std::cout << "  " << cell * 100 / num_cells
                << "%" << std::endl;
    }
#endif
  }  // for(cells)

#ifdef GLENN_DEBUG
  if (chem_->verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "  Chemistry_PK::advance() : "
              << "max iterations - " << max_iterations << " " << "  cell id: " << imax << std::endl;
    std::cout << "  Chemistry_PK::advance() : "
              << "min iterations - " << min_iterations << " " << "  cell id: " << imin << std::endl;
    std::cout << "  Chemistry_PK::advance() : "
              << "ave iterations - " << static_cast<float>(ave_iterations) / num_cells << std::endl;
  }
#endif

  if (chem_->verbosity() == kDebugChemistryProcessKernel) {
    // dumping the values of the final cell. not very helpful by itself,
    // but can be move up into the loops....
    chem_->DisplayTotalColumnHeaders();
    chem_->DisplayTotalColumns(current_time_, beaker_components_);
  }
}  // end advance()


// the MPC will call this function to signal to the
// process kernel that it has accepted the
// state update, thus, the PK should update
// possible auxilary state variables here
void Chemistry_PK::commit_state(Teuchos::RCP<Chemistry_State> chem_state,
                                const double& delta_time) {
  if (chem_->verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "  Chemistry_PK::commit_state() : "
              << "Committing internal state." << std::endl;
  }

  saved_time_ += delta_time;

  if (verbosity() >= kDebugNever) {
    chem_->Speciate(beaker_components_, beaker_parameters_);
    chem_->DisplayResults();
    chem_->DisplayTotalColumnHeaders();
    chem_->DisplayTotalColumns(saved_time_, beaker_components_);
  }
}  // end commit_state()



Teuchos::RCP<Epetra_MultiVector> Chemistry_PK::get_extra_chemistry_output_data() {
  int num_cells = chemistry_state_->porosity()->MyLength();

  for (int cell = 0; cell < num_cells; cell++) {
    // populate aux_data_ by copying from the appropriate internal storage
    // for now, assume we are just looking at free ion conc of primaries
    for (unsigned int i = 0; i < aux_names_.size(); i++) {
      double* cell_aux_data = (*aux_data_)[i];
      double* cell_free_ion = (*chemistry_state_->free_ion_species())[aux_index_.at(i)];
      cell_aux_data[cell] = cell_free_ion[cell];
      if (aux_names_.at(i) == "pH") {
        cell_aux_data[cell] = -std::log10(cell_aux_data[cell]);
      }
    }
  }  // for(cells)

  // return the multi vector
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
  chem_->GetPrimaryNames(names);
}  // end set_component_names()

}  // namespace chemistry
}  // namespace amanzi
