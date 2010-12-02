/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <algorithm>

#include "Chemistry_PK.hpp"
#include "ChemistryException.hpp"
#include "Epetra_MultiVector.h"
#include "SimpleThermoDatabase.hpp"
#include "Beaker.hpp"
#include "Verbosity.hpp"


/*******************************************************************************
 **
 **  Purpose: Trilinos based process kernel for chemistry
 **
 **  Notes:
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
 **  TODO: add method RevertToSavedState() which copies (not swaps)
 **  the saved state back into the current state.
 **
 ******************************************************************************/

Chemistry_PK::Chemistry_PK(Teuchos::ParameterList &param_list,
                           Teuchos::RCP<Chemistry_State> chem_state)
    : status_(ChemistryException::kOkay),
      verbosity_(kSilent),
      max_time_step_(9.9e9),
      chemistry_state_(chem_state),
      parameter_list_(param_list),
      chem_(NULL),
      current_time_(0.0),
      saved_time_(0.0),
      number_aqueous_components_(0),
      number_free_ion_(0),
      number_total_sorbed_(0),
      number_minerals_(0),
      number_ion_exchange_sites_(0),
      number_sorption_sites_(0),
      using_sorption_(false),
      have_free_ion_guess_(false)
{
  // TODO: is there any reason to keep any of this here, or should it
  // all be moved down into InitializeChemistry()...?

  // determine the format of the database file
  std::string database_format =
      parameter_list_.get<std::string>("Thermodynamic Database Format",
                                       "simple");

  // create the appropriate chemistry object
  if (database_format == "simple") {
    chem_ = new SimpleThermoDatabase();
  } else if (database_format == "PFloTran Parsed Database") {
    // hammond stuff here....
  } else {
    // invalid database format, helpful error message and throw an error.
  }

  set_verbosity(static_cast<Verbosity>(parameter_list_.get<int>("Verbosity", 0)));

  //set_verbosity(kDebugChemistryProcessKernel);

  XMLParameters();

  InitializeInternalStorage(&current_state_);
  InitializeInternalStorage(&saved_state_);

  // get initial conditions for minerals etc
  LocalInitialConditions();

  // TODO: finish setting up internal storage for saving previous
  // solutions

}  // end Chemistry_PK()

Chemistry_PK::~Chemistry_PK()
{
  delete chem_;
}  // end ~Chemistry_PK()

// Note: bja: I don't think we want InitializeChemsistry to get a pointer to
// total_component_conc here. We alread have a pointer to the state
// object through chemistry_state_. We should just use that...?
void Chemistry_PK::InitializeChemistry(void)
{
  SizeBeakerComponents();

  // copy the first cell data into the beaker storage for
  // initialization purposes
  int cell=0;
  CopyStateToBeakerParameters(cell);
  CopyCellToBeakerComponents(cell, chemistry_state_->get_total_component_concentration());

  // finish setting up the chemistry object
  try {
    chem_->verbosity(verbosity());
    chem_->Setup(beaker_components_, beaker_parameters_);
    if (verbosity() > kTerse) {
      chem_->Display();
    }

    // solve for initial free-ion concentrations
    chem_->Speciate(beaker_components_, beaker_parameters_);
    if (verbosity() > kTerse) {
      std::cout << "\nTest solution of initial conditions in cell 0:"
                << std::endl;
      chem_->DisplayResults();
    }
  }
  catch (ChemistryException& geochem_error) {
    // TODO: any errer in the constructor is probably fatal. should be
    // rethrown rather than set_status()?
    std::cout << geochem_error.what() << std::endl;
    set_status(geochem_error.error_status());
  }

  CopyBeakerComponentsToCell(cell);

  // now take care of the remainder
  int num_cells = chemistry_state_->get_porosity()->MyLength();
  for (int icell = 1; icell < num_cells; icell++) {
    CopyStateToBeakerParameters(icell);
    CopyCellToBeakerComponents(icell, chemistry_state_->get_total_component_concentration());

    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "Reacting in cell " << icell << std::endl;
    }

    try {
    // solve for initial free-ion concentrations
      chem_->Speciate(beaker_components_, beaker_parameters_);
    }
    catch (ChemistryException& geochem_error) {
      std::cout << geochem_error.what() << std::endl;
      set_status(geochem_error.error_status());
    }
    // if successful copy back
    CopyBeakerComponentsToCell(icell);
  }

} // end InitializeChemistry()

/*******************************************************************************
 **
 **  initialization helper functions
 **
 ******************************************************************************/
void Chemistry_PK::XMLParameters(void)
{
  // extract parameters from the xml list and set in the parameters
  // structure

  if (verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "  Chemistry_PK::XMLParameters()" << std::endl;
  }

  beaker_parameters_ = chem_->GetDefaultParameters();

  // thermo file name
  beaker_parameters_.thermo_database_file =
      parameter_list_.get<std::string>("Thermodynamic Database File",
                                       "dummy.dbs");

  // activity model
  beaker_parameters_.activity_model_name =
      parameter_list_.get<std::string>("Activity Model", "unit");

  // solver parameters here....
  beaker_parameters_.tolerance =
      parameter_list_.get<double>("Tolerance", 1.0e-12);

  // TODO: using <unsigned int> in the parameter list doesn't work...?
  beaker_parameters_.max_iterations =
      static_cast<unsigned int>(parameter_list_.get<int>("Maximum Newton Iterations", 200));

  // other local PK flags
  set_max_time_step(parameter_list_.get<double>("Max Time Step (s)", 9.9e9));

  std::string have_sorption = parameter_list_.get<string>("Using sorption", "no");
  if (have_sorption == "yes") {
    set_using_sorption(true);
  }
  
  // must always have free ions this, the flag only controls if we look in the xml file
  // for an initial guess
  std::string free_ion_guess = parameter_list_.get<string>("Free ion concentrations provided", "no");
  if (free_ion_guess == "yes") {
    set_have_free_ion_guess(true);
  }


}  // end XMLParameters()

// Note: total aqueous components, free ions, total sorbed must all be
// the same size. Each array has it's own number_XXX variable to keep
// life simple. When you want to access total_sorbed, use
// number_total_sorbed(), instead of having to remember which
// variables use number_aqueous_components and which have their own
// variable.
void Chemistry_PK::InitializeInternalStorage(InternalStorage* storage)
{
  if (verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "  Chemistry_PK::InitializeInternalStorage()" << std::endl;
  }
  // physical parameters we don't change, just point to state object
  storage->porosity = chemistry_state_->get_porosity();
  storage->water_saturation = chemistry_state_->get_water_saturation();
  storage->water_density = chemistry_state_->get_water_density();
  storage->volume = chemistry_state_->get_volume();

  // things we need a local copy of and know about because of the
  // state vector
  set_number_aqueous_components(
      chemistry_state_->get_total_component_concentration()->NumVectors());
  storage->aqueous_components = Teuchos::rcp( new Epetra_MultiVector(
      chemistry_state_->get_mesh_maps()->cell_map(false),
      number_aqueous_components() ));
  
  set_number_free_ion(number_aqueous_components());
  storage->free_ion_species = Teuchos::rcp( new Epetra_MultiVector(
            chemistry_state_->get_mesh_maps()->cell_map(false),
            number_free_ion() ));

  // don't know yet if we have these... need to look in the xml
  // file... the main state object really should be reading these
  // in....
  if (parameter_list_.isSublist("Initial Conditions")) {
    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "  Chemistry_PK::InitializeInternalStorage() : "
                << "found initial conditions block in xml data." << std::endl;
    }
    Teuchos::ParameterList initial_conditions =
        parameter_list_.sublist("Initial Conditions");

    set_number_minerals(
        initial_conditions.get<int>("Number of minerals", 0));

    set_number_ion_exchange_sites(
        initial_conditions.get<int>("Number of ion exchange sites", 0));

    set_number_sorption_sites(
        initial_conditions.get<int>("Number of sorption sites", 0));

  }

  if (number_minerals() > 0) {
    storage->minerals = Teuchos::rcp( new Epetra_MultiVector(
        chemistry_state_->get_mesh_maps()->cell_map(false),
        number_minerals() ));
  }

  if (number_ion_exchange_sites()) {
    storage->ion_exchange_sites = Teuchos::rcp( new Epetra_MultiVector(
        chemistry_state_->get_mesh_maps()->cell_map(false),
        number_ion_exchange_sites() ));
  }


  if (using_sorption() == true) {
    set_number_total_sorbed(number_aqueous_components());

    storage->total_sorbed = Teuchos::rcp( new Epetra_MultiVector(
        chemistry_state_->get_mesh_maps()->cell_map(false),
        number_total_sorbed() ));

    if (number_sorption_sites() > 0) {
      storage->sorption_sites = Teuchos::rcp( new Epetra_MultiVector(
              chemistry_state_->get_mesh_maps()->cell_map(false),
              number_sorption_sites() ));
    }
  }
}  // end InitializeInternalStorage()

void Chemistry_PK::SwapCurrentAndSavedStorage(void)
{
  InternalStorage temp;

  std::swap(current_state_.porosity, saved_state_.porosity);

  std::swap(current_state_.water_saturation, 
            saved_state_.water_saturation);

  std::swap(current_state_.water_density, 
            saved_state_.water_density);

  std::swap(current_state_.volume, saved_state_.volume);

  std::swap(current_state_.aqueous_components, 
            saved_state_.aqueous_components);

  std::swap(current_state_.free_ion_species, 
            saved_state_.free_ion_species);

  std::swap(current_state_.minerals, saved_state_.minerals);

  std::swap(current_state_.ion_exchange_sites, 
            saved_state_.ion_exchange_sites);

  std::swap(current_state_.sorption_sites, 
            saved_state_.sorption_sites);

  std::swap(current_state_.total_sorbed, 
            saved_state_.total_sorbed);

  //  end SwapCurrentAndSavedStorage()
}

void Chemistry_PK::LocalInitialConditions(void)
{
  if (verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "  Chemistry_PK::LocalInitialConditions()" << std::endl;
  }
  // add the initial conditions for chemistry specific components:
  // minerals, ion exchange, surface complexation, etc....

  // chem hasn't read the input file yet, so we don't know how many
  // mineral components there are yet... it needs to be in the xml
  // input data....

  // TODO: a lot of duplicate code between here and MPC/State
  // for pulling in initial conditions based on mesh block....

  if (verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "    Looking for initial conditions xml data... ";
  }

  if (parameter_list_.isSublist("Initial Conditions")) {
    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "found." << std::endl;
    }

    Teuchos::ParameterList initial_conditions =
        parameter_list_.sublist("Initial Conditions");

    //
    // loop over mesh blocks
    //

    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "        Expected number of mineral initial conditions: "
                << number_minerals() << std::endl;
      std::cout << "        Expected number of ion exchange sites: "
                << number_ion_exchange_sites() << std::endl;
      std::cout << "        Expected number of sorption sites: "
                << number_sorption_sites() << std::endl;
      if (have_free_ion_guess()) {
        std::cout << "        Expected number of free ion guesses: "
                  << number_free_ion() << std::endl;
      }
    }

    int number_mesh_blocks =
        initial_conditions.get<int>("Number of mesh blocks");

    // TODO: need some sanity check to verify that the number of mesh
    // blocks agrees with the mesh....

    // mesh block count starts at 1!?
    for (int block = 1; block <= number_mesh_blocks; block++) {
      std::stringstream block_name;
      block_name << "Mesh block " << block;
      Teuchos::ParameterList mesh_block_list =
          initial_conditions.sublist(block_name.str());

      int mesh_block_ID = mesh_block_list.get<int>("Mesh block ID");
      if (!chemistry_state_->get_mesh_maps()->valid_set_id(mesh_block_ID,
                                                           Mesh_data::CELL)) {
        // there is an inconsistency in the xml input file...
        throw std::exception();
      }

      //
      // If we have free ion species concentrations in the input file
      // we need to use them (must be of size number_aqueous_components)
      //
      if (have_free_ion_guess()) {
        std::string type("Free Ion Species");
        std::string keyword("Free Ion Species");
        int number_to_find = number_free_ion();
        ExtractInitialCondition(type, keyword, number_to_find, block,
                                mesh_block_list, mesh_block_ID,
                                current_state_.free_ion_species);
        ExtractInitialCondition(type, keyword, number_to_find, block,
                                  mesh_block_list, mesh_block_ID,
                                  saved_state_.free_ion_species);
      } else {
        // need to manually add an initial condition
        std::vector<double> values(number_free_ion(), 1.0e-9);
        set_const_values_for_block(values, number_free_ion(),
                                   current_state_.free_ion_species, mesh_block_ID); 
        set_const_values_for_block(values, number_free_ion(),
                                   saved_state_.free_ion_species, mesh_block_ID); 
      }

      //
      // look for minerals in this mesh block
      //
      if (number_minerals() > 0) {
        std::string type("Minerals");
        std::string keyword("Mineral");
        int number_to_find = number_minerals();
        ExtractInitialCondition(type, keyword, number_to_find, block,
                                mesh_block_list, mesh_block_ID,
                                current_state_.minerals);
        ExtractInitialCondition(type, keyword, number_to_find, block,
                                mesh_block_list, mesh_block_ID,
                                saved_state_.minerals);

      }

      //
      // look for ion exchange sites in this mesh block
      //
      if (number_ion_exchange_sites() > 0) {
        std::string type("Ion Exchange Sites");
        std::string keyword("Ion Exchange Site");
        int number_to_find = number_ion_exchange_sites();
        ExtractInitialCondition(type, keyword, number_to_find, block,
                                mesh_block_list, mesh_block_ID,
                                current_state_.ion_exchange_sites);
        ExtractInitialCondition(type, keyword, number_to_find, block,
                                mesh_block_list, mesh_block_ID,
                                saved_state_.ion_exchange_sites);
      }

      //
      // look for sorption/surface complexation sites
      //
      if (number_sorption_sites() > 0) {
        std::string type("Sorption Sites");
        std::string keyword("Sorption Site");
        int number_to_find = number_sorption_sites();
        ExtractInitialCondition(type, keyword, number_to_find, block,
                                mesh_block_list, mesh_block_ID,
                                current_state_.sorption_sites);
        ExtractInitialCondition(type, keyword, number_to_find, block,
                                mesh_block_list, mesh_block_ID,
                                saved_state_.sorption_sites);
      }
    }  // for (mesh blocks)
  }  // if(initial conditions)

}  // end LocalInitialConditions()

void Chemistry_PK::ExtractInitialCondition(
    const std::string& type,
    const std::string& keyword,
    const int number_to_find,
    const int block,
    const Teuchos::ParameterList& mesh_block_list,
    const int mesh_block_ID,
    Teuchos::RCP<Epetra_MultiVector> data)
{
  if (verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "      Looking for \'" << type
              << "\' initial conditions in block "
              << block << "... ";
  }
  if (mesh_block_list.isSublist(type)) {
    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "found." << std::endl;
    }

    Teuchos::ParameterList ic_list = mesh_block_list.sublist(type);

    std::vector<double> found_values(number_to_find);

    for (int n = 0; n < number_to_find; n++) {
      std::stringstream ic_name;
      ic_name << keyword << " " << n;
      found_values[n] = ic_list.get<double>(ic_name.str());
      if (verbosity() == kDebugChemistryProcessKernel) {
        std::cout << "        Added initial condition for: " << ic_name.str()
                  << "   ic: " << found_values[n] << std::endl;
      }
    }  // for(n)

    set_const_values_for_block(found_values, number_to_find,
                               data, mesh_block_ID);
  }  // if(found type)
  else {
    std::cout << "none." << std::endl;
  }
}  // end ExtractInitialConditions


void Chemistry_PK::set_const_values_for_block(
    const std::vector<double>& values,
    const int num_values,
    Teuchos::RCP<Epetra_MultiVector>& multi_vec,
    const int mesh_block_ID)
{
  for (int n=0; n < num_values; n++) {
    set_cell_value_in_mesh_block(values.at(n), *(*multi_vec)(n),
                                 mesh_block_ID);
  }
}  // end set_const_values_for_block()


void Chemistry_PK::set_cell_value_in_mesh_block(const double value,
                                                Epetra_Vector& vec,
                                                const int mesh_block_id)
{
  if (!chemistry_state_->get_mesh_maps()->valid_set_id(mesh_block_id,
                                                       Mesh_data::CELL)) {
    throw std::exception();
  }

  unsigned int mesh_block_size =
      chemistry_state_->get_mesh_maps()->get_set_size(mesh_block_id,
                                                      Mesh_data::CELL,
                                                      OWNED);

  std::vector<unsigned int> cell_ids(mesh_block_size);

  chemistry_state_->get_mesh_maps()->get_set(mesh_block_id,
                                             Mesh_data::CELL, OWNED,
                                             cell_ids.begin(), cell_ids.end());

  for( std::vector<unsigned int>::iterator c = cell_ids.begin();
       c != cell_ids.end();  c++) {
    vec[*c] = value;
  }

}  // end set_cell_value_in_mesh_block()


void Chemistry_PK::SizeBeakerComponents(void)
{
  // initialize the beaker component data structure
  beaker_components_.total.clear();
  beaker_components_.free_ion.clear();
  beaker_components_.minerals.clear();
  beaker_components_.ion_exchange_sites.clear();
  beaker_components_.total_sorbed.clear();

  beaker_components_.total.resize(number_aqueous_components(), 0.0);
  beaker_components_.minerals.resize(number_minerals(), 0.0);
  beaker_components_.ion_exchange_sites.resize(number_ion_exchange_sites(), 0.0);
  if (using_sorption()) {
    beaker_components_.total_sorbed.resize(number_total_sorbed(), 0.0);
    // sorption sites would go here....
  }
  // free_ion needs a non-zero default value if no initial guess provide in input
  beaker_components_.free_ion.resize(number_aqueous_components(), 1.0e-9);
}  // end SizeBeakerComponents()


void Chemistry_PK::CopyCellToBeakerComponents(
    const int cell_id,
    Teuchos::RCP<const Epetra_MultiVector> aqueous_components)
{
  // copy component data from the cell arrays into beaker component
  // structure

  // Note: total uses the passed in value from transport, not the
  // current storage value
  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_components = (*aqueous_components)[c];
    beaker_components_.total[c] = cell_components[cell_id];
  }

  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_free_ion = (*current_state_.free_ion_species)[c];
    beaker_components_.free_ion[c] = cell_free_ion[cell_id];
  }

  for (unsigned int m = 0; m < number_minerals(); m++) {
    double* cell_minerals = (*current_state_.minerals)[m];
    beaker_components_.minerals[m] = cell_minerals[cell_id];
  }

  for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) {
    double* cell_ion_exchange_sites = (*current_state_.ion_exchange_sites)[i];
    beaker_components_.ion_exchange_sites[i] = cell_ion_exchange_sites[cell_id];
  }

  if (using_sorption()) {
    for (unsigned int c = 0; c < number_aqueous_components(); c++) {
      double* cell_total_sorbed = (*current_state_.total_sorbed)[c];
      beaker_components_.total_sorbed[c] = cell_total_sorbed[cell_id];
    }
//     if (number_sorption_sites() > 0) {
//       for (unsigned int i = 0; i < number_sorption_sites(); i++) {
//         double* cell_sorption_sites = (*current_state_.sorption_sites)[i];
//         beaker_components_.sorption_sites[i] = cell_sorption_sites[cell_id];
//       }
//     }
  }  // end if(using_sorption)

}  // end CopyCellToBeakerComponents()


void Chemistry_PK::CopyBeakerComponentsToCell(const int cell_id)
{
  // copy data from the beaker back into the state arrays

  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_components = (*current_state_.aqueous_components)[c];
    cell_components[cell_id] = beaker_components_.total[c];
  }

  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_free_ion = (*current_state_.free_ion_species)[c];
    cell_free_ion[cell_id] = beaker_components_.free_ion[c];
  }

  for (unsigned int m = 0; m < number_minerals(); m++) {
    double* cell_minerals = (*current_state_.minerals)[m];
    cell_minerals[cell_id] = beaker_components_.minerals[m];
  }

  for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) {
    double* cell_ion_exchange_sites = (*current_state_.ion_exchange_sites)[i];
    cell_ion_exchange_sites[cell_id] = beaker_components_.ion_exchange_sites[i];
  }

  if (using_sorption()) {
    for (unsigned int c = 0; c < number_aqueous_components(); c++) {
      double* cell_total_sorbed = (*current_state_.total_sorbed)[c];
      cell_total_sorbed[cell_id] = beaker_components_.total_sorbed[c];
    }
//     if (number_sorption_sites() > 0) {
//       for (unsigned int i = 0; i < number_sorption_sites(); i++) {
//         double* cell_sorption_sites = (*current_state_.sorption_sites)[i];
//         cell_sorption_sites[cell_id] = beaker_components_.total_sorbed[i];
//       }
//     }
  }  // end if(using_sorption)

}  // end CopyBeakerComponentsToCell()


void Chemistry_PK::CopyStateToBeakerParameters(const int cell_id)
{
  // copy data from state arrays into the beaker parameters
  beaker_parameters_.water_density = (*current_state_.water_density)[cell_id];
  beaker_parameters_.porosity = (*current_state_.porosity)[cell_id];
  beaker_parameters_.saturation = (*current_state_.water_saturation)[cell_id];
  beaker_parameters_.volume = (*current_state_.volume)[cell_id];
}  // end CopyStateToBeakerParameters()



/*******************************************************************************
 **
 **  MPC interface functions
 **
 ******************************************************************************/
Teuchos::RCP<Epetra_MultiVector> Chemistry_PK::get_total_component_concentration(void) const
{
  return current_state_.aqueous_components;
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
    Teuchos::RCP<const Epetra_MultiVector> total_component_concentration_star)
{
  if (chem_->verbosity() == kDebugChemistryProcessKernel) {
    cout << "  Chemistry_PK::advance() : "
         << "advancing the chemistry process model..." << endl;
  }

  current_time_ = saved_time_ + delta_time;

  // shorter name for the state that came out of transport
  Teuchos::RCP<const Epetra_MultiVector> tcc_star =
      total_component_concentration_star;


  // TODO: use size of the porosity vector as indicator of size for
  // now... should get data from the mesh...?
  int num_cells = chemistry_state_->get_porosity()->MyLength();

  for (int cell = 0; cell < num_cells; cell++) {
    // copy state data into the beaker data structures
    CopyCellToBeakerComponents(cell, tcc_star);
    CopyStateToBeakerParameters(cell);

    try {
      // chemistry computations for this cell
      chem_->ReactionStep(&beaker_components_, beaker_parameters_, delta_time);
    }
    catch (ChemistryException& geochem_error) {
      std::cout << "ERROR: Chemistry_PR::advance() "
                << "cell[" << cell << "]: " << std::endl;
      std::cout << geochem_error.what();
      chem_->DisplayTotalColumnHeaders();
      chem_->DisplayTotalColumns(current_time_, beaker_components_.total);
      std::cout << std::endl;
      set_status(geochem_error.error_status());
    }

    // update this cell's data in the arrays
    CopyBeakerComponentsToCell(cell);

    // TODO: was porosity etc changed? copy someplace

  }  // for(cells)

  if (chem_->verbosity() == kDebugChemistryProcessKernel) {
    // dumping the values of the final cell. not very helpful by itself,
    // but can be move up into the loops....
    chem_->DisplayTotalColumnHeaders();
    chem_->DisplayTotalColumns(current_time_, beaker_components_.total);
  }

}  // end advance()


// the MPC will call this function to signal to the
// process kernel that it has accepted the
// state update, thus, the PK should update
// possible auxilary state variables here
void Chemistry_PK::commit_state(Teuchos::RCP<Chemistry_State> chem_state,
                                const double& delta_time)
{
  if (chem_->verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "  Chemistry_PK::commit_state() : "
              << "Committing internal state." << std::endl;
  }

  saved_time_ += delta_time;

  SwapCurrentAndSavedStorage();

  if (verbosity() >= kTerse) {
    //chem_->Speciate(beaker_components_, beaker_parameters_);
    //chem_->DisplayResults();
    //chem_->DisplayTotalColumnHeaders();
    //chem_->DisplayTotalColumns(saved_time_, beaker_components_.total);
  }
}  // end commit_state()



Teuchos::RCP<Epetra_MultiVector> Chemistry_PK::get_extra_chemistry_output_data()
{
  // get the number of extra output variables, for now I'm using 1

  int nvars = 1;
 
  // create the Epetra_MultiVector that will hold the data

  Teuchos::RCP<Epetra_MultiVector> aux_data = 
    Teuchos::rcp(new Epetra_MultiVector(chemistry_state_->get_mesh_maps()->cell_map(false), nvars));

  // populate aux_data, for now I am just initializing to 0, vector by vector
  for (int i=0; i<nvars; i++) 
    {
      // the i-th Epetra_Vecotr in the aux_data Epetre_MultiVector
      // is *(*aux_data)(i)

      (*aux_data)(i)->PutScalar(0.0);
    }

  // return the multi vector
  return aux_data;
}


void Chemistry_PK::set_chemistry_output_names(std::vector<string> &names)
{
  // get the number of auxillary variables, for now, just one

  int nvars = 1;

  // resize the name vector
  names.resize(nvars);
  
  // populate the name vactor, for now there's only one
  for (int i=0; i<nvars; i++)
    {
      names[i] = "auxillary variable";
    }
}


void Chemistry_PK::set_component_names(std::vector<string> &names)
{
  // get the number of components
  int ncomps = chemistry_state_->get_total_component_concentration()->NumVectors();

  // resize the names vector 
  names.resize(ncomps);

  // populate the names vector, for testing I'm just adding "-test"

  for (int i=0; i<ncomps; i++)
    {
      names[i] = "test";
    }
}
