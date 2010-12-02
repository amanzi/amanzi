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
      using_sorption_(0),
      have_free_ion_guess_(0)
{
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

  SizeBeakerComponents();

  // copy the first cell data into the beaker storage for
  // initialization purposes
  int cell=0;
  CopyStateToBeakerParameters(cell);
  CopyCellToBeakerComponents(cell, current_state_.aqueous_components);


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

  // loop through every cell and verify that the initial conditions
  // produce a valid solution...? reaction step or speciate...?

  // TODO: finish setting up internal storage for saving previous
  // solutions

}  // end Chemistry_PK()

Chemistry_PK::~Chemistry_PK()
{
  delete chem_;
}  // end ~Chemistry_PK()


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

  set_max_time_step(parameter_list_.get<double>("Max Time Step (s)", 9.9e9));

}  // end XMLParameters()


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

  // don't know yet if we have these... need to look in the xml
  // file... the main state object really should be reading these
  // in....
  //storage->minerals = Teuchos::ENull;
  //storage->ion_exchange_sites = Teuchos::ENull;
  //storage->sorption_sites = Teuchos::ENull;
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

  std::swap(current_state_.free_ion, 
            saved_state_.free_ion);

  std::swap(current_state_.total_sorbed, 
            saved_state_.total_sorbed);

  std::swap(current_state_.minerals, saved_state_.minerals);

  std::swap(current_state_.ion_exchange_sites, 
            saved_state_.ion_exchange_sites);

  std::swap(current_state_.sorption_sites, 
            saved_state_.sorption_sites);

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

  //
  // do we need have sorption?
  //
  set_using_sorption(parameter_list_.get<int>("Using sorption", 0));
  if (using_sorption()) {
    set_number_total_sorbed(number_aqueous_components());
    current_state_.total_sorbed =
        Teuchos::rcp( new Epetra_MultiVector(
            chemistry_state_->get_mesh_maps()->cell_map(false),
            number_total_sorbed() ));
    saved_state_.total_sorbed =
        Teuchos::rcp( new Epetra_MultiVector(
            chemistry_state_->get_mesh_maps()->cell_map(false),
            number_total_sorbed() ));
  }

  if (parameter_list_.isSublist("Initial Conditions")) {
    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "found." << std::endl;
    }

    Teuchos::ParameterList initial_conditions =
        parameter_list_.sublist("Initial Conditions");

    //
    // do we have minerals
    //
    set_number_minerals(
        initial_conditions.get<int>("Number of minerals", 0));

    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "        Expected number of mineral initial conditions: "
                << number_minerals() << std::endl;
    }

    if (number_minerals()) {
      current_state_.minerals =
          Teuchos::rcp( new Epetra_MultiVector(
              chemistry_state_->get_mesh_maps()->cell_map(false),
              number_minerals() ));
      saved_state_.minerals =
          Teuchos::rcp( new Epetra_MultiVector(
              chemistry_state_->get_mesh_maps()->cell_map(false),
              number_minerals() ));
    }

    //
    // do we have a guess for free ions?
    //
    set_have_free_ion_guess(
        initial_conditions.get<int>("Have free ion guess", 0));
    // must always set this, the flag only controls if we look in the xml file
    set_number_free_ion(number_aqueous_components());
    if (number_free_ion()) {
      current_state_.free_ion =
          Teuchos::rcp( new Epetra_MultiVector(
              chemistry_state_->get_mesh_maps()->cell_map(false),
              number_free_ion() ));
      saved_state_.free_ion =
          Teuchos::rcp( new Epetra_MultiVector(
              chemistry_state_->get_mesh_maps()->cell_map(false),
              number_free_ion() ));
    }

    //
    // do we have ion exchange sites
    //

    set_number_ion_exchange_sites(
        initial_conditions.get<int>("Number of ion exchange sites", 0));
    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "        Expected number of ion exchange sites: "
                << number_ion_exchange_sites() << std::endl;
    }

    if (number_ion_exchange_sites()) {
      current_state_.ion_exchange_sites =
          Teuchos::rcp( new Epetra_MultiVector(
              chemistry_state_->get_mesh_maps()->cell_map(false),
              number_ion_exchange_sites() ));
      saved_state_.ion_exchange_sites =
          Teuchos::rcp( new Epetra_MultiVector(
              chemistry_state_->get_mesh_maps()->cell_map(false),
              number_ion_exchange_sites() ));
    }

    //
    // do we have surface complexes
    //

    set_number_sorption_sites(
        initial_conditions.get<int>("Number of sorption sites", 0));

    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "        Expected number of sorption sites: "
                << number_sorption_sites() << std::endl;
    }

    if (number_sorption_sites()) {
      current_state_.sorption_sites =
          Teuchos::rcp( new Epetra_MultiVector(
              chemistry_state_->get_mesh_maps()->cell_map(false),
              number_sorption_sites() ));
      saved_state_.sorption_sites =
          Teuchos::rcp( new Epetra_MultiVector(
              chemistry_state_->get_mesh_maps()->cell_map(false),
              number_sorption_sites() ));
    }

    //
    // loop over mesh blocks
    //

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
      // look for free ions
      //
      if (have_free_ion_guess()) {
        std::string type("Free ions");
        std::string keyword("Component");
        int number_to_find = number_free_ion();
        ExtractInitialCondition(type, keyword, number_to_find, block,
                                mesh_block_list, mesh_block_ID,
                                current_state_.free_ion);
        ExtractInitialCondition(type, keyword, number_to_find, block,
                                mesh_block_list, mesh_block_ID,
                                saved_state_.free_ion);

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

  beaker_components_.total.resize(number_aqueous_components());
  beaker_components_.free_ion.resize(number_aqueous_components());
  beaker_components_.minerals.resize(number_minerals());
  beaker_components_.ion_exchange_sites.resize(number_ion_exchange_sites());
  beaker_components_.total_sorbed.resize(number_total_sorbed());
}  // end SizeBeakerComponents()


void Chemistry_PK::CopyCellToBeakerComponents(
    const int cell_id,
    Teuchos::RCP<const Epetra_MultiVector> aqueous_components)
{
  // copy component data from the cell arrays into beaker component
  // structure

  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_components = (*aqueous_components)[c];
    beaker_components_.total[c] = cell_components[cell_id];
  }

  for (unsigned int f = 0; f < number_free_ion(); f++) {
    double* cell_free_ion = (*current_state_.free_ion)[f];
    beaker_components_.free_ion[f] = cell_free_ion[cell_id];
  }

  for (unsigned int m = 0; m < number_minerals(); m++) {
    double* cell_minerals = (*current_state_.minerals)[m];
    beaker_components_.minerals[m] = cell_minerals[cell_id];
  }

  for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) {
    double* cell_ion_exchange_sites = (*current_state_.ion_exchange_sites)[i];
    beaker_components_.ion_exchange_sites[i] = cell_ion_exchange_sites[cell_id];
  }

  for (unsigned int i = 0; i < number_total_sorbed(); i++) {
    double* cell_total_sorbed = (*current_state_.total_sorbed)[i];
    beaker_components_.total_sorbed[i] = cell_total_sorbed[cell_id];
  }

//   for (unsigned int i = 0; i < number_sorption_sites(); i++) {
//     double* cell_sorption_sites = (*current_state_.sorption_sites)[i];
//     beaker_components_.total_sorbed[i] = cell_sorption_sites[cell_id];
//   }

}  // end CopyCellToBeakerComponents()


void Chemistry_PK::CopyBeakerComponentsToCell(const int cell_id)
{
  // copy data from the beaker back into the state arrays

  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_components = (*current_state_.aqueous_components)[c];
    cell_components[cell_id] = beaker_components_.total[c];
  }

  for (unsigned int m = 0; m < number_minerals(); m++) {
    double* cell_minerals = (*current_state_.minerals)[m];
    cell_minerals[cell_id] = beaker_components_.minerals[m];
  }

  for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) {
    double* cell_ion_exchange_sites = (*current_state_.ion_exchange_sites)[i];
    cell_ion_exchange_sites[cell_id] = beaker_components_.ion_exchange_sites[i];
  }

  for (unsigned int i = 0; i < number_sorption_sites(); i++) {
    double* cell_sorption_sites = (*current_state_.sorption_sites)[i];
    cell_sorption_sites[cell_id] = beaker_components_.total_sorbed[i];
  }

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
