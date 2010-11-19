/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
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
 *******************************************************************************/

Chemistry_PK::Chemistry_PK(Teuchos::ParameterList &param_list,
                           Teuchos::RCP<Chemistry_State> chem_state)
    : status_(ChemistryException::kOkay),
      verbosity_(kSilent),
      chemistry_state_(chem_state),
      parameter_list_(param_list),
      chem_(NULL),
      current_time_(0.0),
      saved_time_(0.0),
      number_aqueous_components_(0),
      number_minerals_(0),
      number_ion_exchange_sites_(0),
      number_sorption_sites_(0)
{
  // determine the format of the database file
  std::string database_format =
      parameter_list_.get<std::string>("Thermodynamic Database Format", "simple");

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

  LocalPhysicalState();

  // assume that the State object read in appropriate initial
  // conditions for the total components, trac #216, comment 7
  Teuchos::RCP<const Epetra_MultiVector> temp = chemistry_state_->get_total_component_concentration();
  aqueous_components_ = Teuchos::rcp(new Epetra_MultiVector(*temp));

  set_number_aqueous_components(aqueous_components_->NumVectors());

  // get initial conditions for minerals etc
  LocalInitialConditions();

  SizeBeakerComponents();

  // copy the first cell data into the beaker storage for
  // initialization purposes
  int cell=0;
  CopyCellToBeakerComponents(cell, aqueous_components_);


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
      std::cout << "\nTesting initial solution in cell 0:" << std::endl;
      chem_->DisplayResults();
    }
  }
  catch (ChemistryException& geochem_error) {
    std::cout << geochem_error.what() << std::endl;
    set_status(geochem_error.error_status());
  }
  // loop through every cell and verify that the initial conditions
  // produce a valid solution...? reaction step or speciate...?

}  // end Chemistry_PK()

Chemistry_PK::~Chemistry_PK()
{
  delete chem_;
}  // end ~Chemistry_PK()


/*******************************************************************************
 **
 **  initialization helper functions
 **
 *******************************************************************************/
void Chemistry_PK::XMLParameters(void)
{
  if (verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "  Chemistry_PK::XMLParameters()" << std::endl;
  }
  // extract and set various parameters
  beaker_parameters_ = chem_->GetDefaultParameters();

  // thermo file name
  beaker_parameters_.thermo_database_file =
      parameter_list_.get<std::string>("Thermodynamic Database File", "dummy.dbs");

  // other input data file name?

  // activity model
  beaker_parameters_.activity_model_name =
      parameter_list_.get<std::string>("Activity Model", "unit");

  // solver parameters here....
  beaker_parameters_.tolerance =
      parameter_list_.get<double>("Tolerance", 1.1e-12);

  // TODO: using <unsigned int> in the parameter list doesn't work...?
  beaker_parameters_.max_iterations =
      static_cast<unsigned int>(parameter_list_.get<int>("Maximum Newton Iterations", 200));

}  // end XMLParameters()

void Chemistry_PK::LocalPhysicalState(void)
{
  if (verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "  Chemistry_PK::LocalPhysicalState()" << std::endl;
  }
  // physical parameters
  current_porosity_ = chemistry_state_->get_porosity();
  beaker_parameters_.porosity = (*current_porosity_)[0];

  current_water_saturation_ = chemistry_state_->get_water_saturation();
  beaker_parameters_.saturation = (*current_water_saturation_)[0];

  current_water_density_ = chemistry_state_->get_water_density();
  beaker_parameters_.water_density = (*current_water_density_)[0];

}  // end LocalPhysicalState()

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

    Teuchos::ParameterList initial_conditions = parameter_list_.sublist("Initial Conditions");

    //
    // do we have minerals, ion exchange sites, surface complexes, etc?
    //
    set_number_minerals(initial_conditions.get<int>("Number of minerals"));
    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "        Expected number of mineral initial conditions: "
                << number_minerals() << std::endl;
    }
    if (number_minerals()) {
      minerals_ = Teuchos::rcp( new Epetra_MultiVector( chemistry_state_->get_mesh_maps()->cell_map(false),
                                                        number_minerals() ));
    }


    set_number_ion_exchange_sites(initial_conditions.get<int>("Number of ion exchange sites", 0));
    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "        Expected number of ion exchange sites: "
                << number_ion_exchange_sites() << std::endl;
    }
    if (number_ion_exchange_sites()) {
      ion_exchange_sites_ =
          Teuchos::rcp( new Epetra_MultiVector( chemistry_state_->get_mesh_maps()->cell_map(false),
                                                number_ion_exchange_sites() ));
    }


    set_number_sorption_sites(initial_conditions.get<int>("Number of sorption sites", 0));
    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "        Expected number of sorption sites: "
                << number_sorption_sites() << std::endl;
    }
    if (number_sorption_sites()) {
      sorption_sites_ =
          Teuchos::rcp( new Epetra_MultiVector( chemistry_state_->get_mesh_maps()->cell_map(false),
                                                number_sorption_sites() ));
    }

    //
    // loop over mesh blocks
    //
    int number_mesh_blocks = initial_conditions.get<int>("Number of mesh blocks");

    // need some sanity check to verify that the number of mesh blocks
    // agrees with the mesh....
    for (int block = 1; block <= number_mesh_blocks; block++) {
      std::stringstream block_name;
      block_name << "Mesh block " << block;
      Teuchos::ParameterList mesh_block_list = initial_conditions.sublist(block_name.str());

      int mesh_block_ID = mesh_block_list.get<int>("Mesh block ID");
      if (!chemistry_state_->get_mesh_maps()->valid_set_id(mesh_block_ID, Mesh_data::CELL)) {
        // there is an inconsistency in the xml input file, die
        throw std::exception();
      }

      //
      // look for minerals in this block
      //
      if (number_minerals() > 0) {
        std::string type("Minerals");
        std::string keyword("Mineral");
        int number_to_find = number_minerals();
        ExtractInitialCondition(type, keyword, number_to_find, block, 
                                mesh_block_list, mesh_block_ID,
                                minerals_);

      }  // if (expect minerals)

      //
      // look for ion exchange sites in this block
      //
      if (number_ion_exchange_sites() > 0) {
        std::string type("Ion Exchange Sites");
        std::string keyword("Ion Exchange Site");
        int number_to_find = number_ion_exchange_sites();
        ExtractInitialCondition(type, keyword, number_to_find, block, 
                                mesh_block_list, mesh_block_ID,
                                ion_exchange_sites_);
      }  // if (expect ion exchange sites)

      //
      // look for sorption/surface complexation sites in this block
      //
      if (number_sorption_sites() > 0) {
        std::string type("Sorption Sites");
        std::string keyword("Sorption Site");
        int number_to_find = number_sorption_sites();
        ExtractInitialCondition(type, keyword, number_to_find, block, 
                                mesh_block_list, mesh_block_ID,
                                sorption_sites_);
      }  // if(expect sorption sites)
    }  // for (mesh blocks)
  }  // if(initial conditions)

}  // end LocalInitialConditions()

void Chemistry_PK::ExtractInitialCondition(const std::string& type,
                                           const std::string& keyword,
                                           const int number_to_find,
                                           const int block,
                                           const Teuchos::ParameterList& mesh_block_list,
                                           const int mesh_block_ID,
                                           Teuchos::RCP<Epetra_MultiVector> data)
{
  if (verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "      Looking for \'" << type << "\' initial conditions in block "
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
    set_const_values_for_block(found_values, number_to_find, data, mesh_block_ID);
  }  // if(found type)
}  // end ExtractInitialConditions


void Chemistry_PK::set_const_values_for_block(const std::vector<double>& values,
                                              const int num_values,
                                              Teuchos::RCP<Epetra_MultiVector>& multi_vec,
                                              const int mesh_block_ID)
{
  for (int n=0; n < num_values; n++) {
    set_cell_value_in_mesh_block(values.at(n), *(*multi_vec)(n), mesh_block_ID);
  }
}

void Chemistry_PK::set_cell_value_in_mesh_block(const double value,
                                                Epetra_Vector& vec,
                                                const int mesh_block_id)
{
  if (!chemistry_state_->get_mesh_maps()->valid_set_id(mesh_block_id, Mesh_data::CELL)) {
    throw std::exception();
  }

  unsigned int mesh_block_size =
      chemistry_state_->get_mesh_maps()->get_set_size(mesh_block_id,
                                                      Mesh_data::CELL,
                                                      OWNED);

  std::vector<unsigned int> cell_ids(mesh_block_size);

  chemistry_state_->get_mesh_maps()->get_set(mesh_block_id, Mesh_data::CELL, OWNED,
                                             cell_ids.begin(), cell_ids.end());

  for( std::vector<unsigned int>::iterator c = cell_ids.begin();
       c != cell_ids.end();  c++) {
    vec[*c] = value;
  }

}

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
  beaker_components_.total_sorbed.resize(number_sorption_sites());
}  // end SizeBeakerComponents()


void Chemistry_PK::CopyCellToBeakerComponents(int cell_id,
                                              Teuchos::RCP<const Epetra_MultiVector> aqueous_components)
{
  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_components = (*aqueous_components)[c];
    beaker_components_.total[c] = cell_components[cell_id];
    beaker_components_.free_ion[c] = 1.0e-9;
  }

  for (unsigned int m = 0; m < number_minerals(); m++) {
    double* cell_minerals = (*minerals_)[m];
    beaker_components_.minerals[m] = cell_minerals[cell_id];
  }

  for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) {
    double* cell_ion_exchange_sites = (*ion_exchange_sites_)[i];
    beaker_components_.ion_exchange_sites[i] = cell_ion_exchange_sites[cell_id];
  }

  for (unsigned int i = 0; i < number_sorption_sites(); i++) {
    double* cell_sorption_sites = (*sorption_sites_)[i];
    beaker_components_.total_sorbed[i] = cell_sorption_sites[cell_id];
  }

}  // end CopyCellToBeakerComponents()

void Chemistry_PK::CopyBeakerComponentsToCell(int cell_id)
{
  for (unsigned int c = 0; c < number_aqueous_components(); c++) {
    double* cell_components = (*aqueous_components_)[c];
    cell_components[cell_id] = beaker_components_.total[c];
  }

  for (unsigned int m = 0; m < number_minerals(); m++) {
    double* cell_minerals = (*minerals_)[m];
    cell_minerals[cell_id] = beaker_components_.minerals[m];
  }

  for (unsigned int i = 0; i < number_ion_exchange_sites(); i++) {
    double* cell_ion_exchange_sites = (*ion_exchange_sites_)[i];
    cell_ion_exchange_sites[cell_id] = beaker_components_.ion_exchange_sites[i];
  }

  for (unsigned int i = 0; i < number_sorption_sites(); i++) {
    double* cell_sorption_sites = (*sorption_sites_)[i];
    cell_sorption_sites[cell_id] = beaker_components_.total_sorbed[i];
  }

}  // end CopyBeakerComponentsToCell()

/*******************************************************************************
 **
 **  MPC interface functions
 **
 *******************************************************************************/
Teuchos::RCP<Epetra_MultiVector> Chemistry_PK::get_total_component_concentration(void) const
{
  return aqueous_components_;
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
    cout << "  Chemistry_PK::advance() : advancing the chemistry process model..." << endl;
  }

  current_time_ = saved_time_ + delta_time;
  Teuchos::RCP<const Epetra_MultiVector> tcc_star = total_component_concentration_star;
  // loop over all cells
  // use size of the porosity vector as indicator of size for now....
  int num_cells = chemistry_state_->get_porosity()->MyLength();
  for (int cell = 0; cell < num_cells; cell++) {
    CopyCellToBeakerComponents(cell, tcc_star);

    // chemistry computations for this cell
    try {
      chem_->ReactionStep(&beaker_components_, beaker_parameters_, delta_time);
    }
    catch (ChemistryException& geochem_error) {
      std::cout << geochem_error.what() << std::endl;
      set_status(geochem_error.error_status());
    }
    // update chemistry
    CopyBeakerComponentsToCell(cell);

  }  // for(cells)


  if (chem_->verbosity() == kDebugChemistryProcessKernel) {
    // dumping the values of the final cell. not very helpful by itself,
    // but can be move up into the loops....
    chem_->DisplayTotalColumns(current_time_, beaker_components_.total);
  }

}  // end advance()


void Chemistry_PK::commit_state(Teuchos::RCP<Chemistry_State> chem_state,
                                const double& delta_time)
{
  if (chem_->verbosity() == kDebugChemistryProcessKernel) {
    cout << "  Chemistry_PK::commit_state() : Committing internal state." << endl;
  }
  // the MPC will call this function to signal to the
  // process kernel that it has accepted the
  // state update, thus, the PK should update
  // possible auxilary state variables here
  saved_time_ += delta_time;
}  // end commit_state()
