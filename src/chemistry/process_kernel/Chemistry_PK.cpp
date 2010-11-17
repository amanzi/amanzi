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
      num_aqueous_components_(0)
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

  // initialize the beaker component data structure
  beaker_components_.free_ion.clear();
  beaker_components_.minerals.clear();
  beaker_components_.ion_exchange_sites.clear();
  beaker_components_.total.clear();
  beaker_components_.total_sorbed.clear();

  // assume that the State object read in appropriate initial
  // conditions for the total components, trac #216, comment 7
  Teuchos::RCP<const Epetra_MultiVector> aqueous_components = 
      chemistry_state_->get_total_component_concentration();
  num_aqueous_components(aqueous_components->NumVectors());

  // copy state into internal storage...?
  
  // copy the first cell data into the beaker storage for
  // initialization purposes
  for (unsigned int c = 0; c < num_aqueous_components(); c++) {
    double* component = (*aqueous_components)[c];
    beaker_components_.total.push_back(component[0]);
  }

  LocalInitialConditions();

  if (beaker_components_.free_ion.size() != beaker_components_.total.size()) {
    beaker_components_.free_ion.resize(beaker_components_.total.size(), 1.0e-9);
  }

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

  // TODO: we are assuming constant local initial conditions over the
  // entire domain. Need to look at the MPC and apply initial
  // conditions for each mesh block group

  if (verbosity() == kDebugChemistryProcessKernel) {
    std::cout << "    Looking for initial conditions xml data... ";
  }
  if (parameter_list_.isSublist("Initial Conditions")) {
    if (verbosity() == kDebugChemistryProcessKernel) {
      std::cout << "found." << std::endl;
    }
    Teuchos::ParameterList initial_conditions = parameter_list_.sublist("Initial Conditions");
    if (initial_conditions.isSublist("Minerals")) {
      if (verbosity() == kDebugChemistryProcessKernel) {
        std::cout << "      Found mineral initial conditions." << std::endl;        
      }
      Teuchos::ParameterList initial_conditions_minerals = initial_conditions.sublist("Minerals");
      int num_minerals = initial_conditions_minerals.get<int>("Number of minerals");
      if (verbosity() == kDebugChemistryProcessKernel) {
        std::cout << "        Expected number of mineral initial conditions: " << num_minerals << std::endl;        
      }      
      beaker_components_.minerals.resize(num_minerals);
      for (int m = 1; m <= num_minerals; m++) {
        // indexing must start from one... see trac# 216, comment 7
        std::stringstream mineral_ic_name;
        mineral_ic_name << "Mineral " << m;
        double mineral_ic = initial_conditions_minerals.get<double>(mineral_ic_name.str(), 0.0);
        // TODO: this needs to be put into an array, size(num_cells)
        beaker_components_.minerals[m] = mineral_ic;
        if (verbosity() == kDebugChemistryProcessKernel) {
          std::cout << "        Added initial condition for: " << mineral_ic_name.str() 
                    << "   ic: " << mineral_ic << std::endl;
        }
      }  // for(m in minerals)
    }  // if(mineral initial conditions)

    if (initial_conditions.isSublist("Ion Exchange")) {
      if (verbosity() == kDebugChemistryProcessKernel) {
        std::cout << "    Chemistry_PK:: found ion exchange initial conditions" << std::endl;        
      }
      std::cout << "    Chemistry_PK:: ion exchange initial conditions not implemented yet." << std::endl;
    }  // if(ion exchange initial conditions)

    if (initial_conditions.isSublist("Surface Complexation")) {
      if (verbosity() == kDebugChemistryProcessKernel) {
        std::cout << "    Chemistry_PK:: found surface complexation initial conditions" << std::endl;        
      }
      std::cout << "    Chemistry_PK:: surface complexation initial conditions not implemented yet." << std::endl;
    }  // if(surface complexation initial conditions)
    
  }  // if(initial conditions)
  
}  // end LocalInitialConditions()


/*******************************************************************************
 **
 **  MPC interface functions
 **
 *******************************************************************************/

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
    Teuchos::RCP<Epetra_MultiVector> total_component_concentration_star)
{
  if (chem_->verbosity() == kDebugChemistryProcessKernel) {
    cout << "  Chemistry_PK::advance() : advancing the chemistry process model..." << endl;
  }

  current_time_ = saved_time_ + delta_time;
  Teuchos::RCP<const Epetra_MultiVector> state_tcc = chemistry_state_->get_total_component_concentration();
  // loop over all cells
  // use size of the porosity vector as indicator of size for now....
  int num_cells = chemistry_state_->get_porosity()->MyLength();
  for (int cell = 0; cell < num_cells; cell++) {

    // TODO: copy the components for this cell out of tcc or tcc_star?
    // tcc* is empty at the moment
    for (unsigned int c = 0; c < num_aqueous_components(); c++) {
      double* component = (*state_tcc)[c];
      beaker_components_.total[c] = component[cell];
    }

    // copy the minerals, surface complexation sites, ion exchange sites, etc
    


    // chemistry computations for this cell
    chem_->ReactionStep(&beaker_components_, beaker_parameters_, delta_time);

    // TODO: extract the results and save them into...? do we really want to
    // over write tcc* if it has the transport solution in it?
    //std::cout << "num vectors in size tcc* = " << total_component_concentration_star->NumVectors() << std::endl; 
    for (unsigned int c = 0; c < num_aqueous_components(); c++) {
      double* component = (*total_component_concentration_star)[c];
      component[cell] = beaker_components_.total[c];
    }

    // update the minerals, surface complexation sites, ion exchange sites, etc



  }
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
