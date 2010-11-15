/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "Chemistry_PK.hpp"
#include "Epetra_MultiVector.h"
#include "SimpleThermoDatabase.hpp"
#include "Beaker.hpp"
#include "Verbosity.hpp"


/*
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
**
*/

Chemistry_PK::Chemistry_PK(Teuchos::ParameterList &param_list,
                           Teuchos::RCP<Chemistry_State> chem_state)
    : status_(kChemistryOK),
      verbosity_(kSilent),
      chemistry_state_(chem_state),
      parameter_list_(param_list),
      chem_(NULL),
      current_time_(0.0),
      saved_time_(0.0),
      num_aqueous_components_(0)
{
  // determine the format of the database file
  std::string database_format = parameter_list_.get<std::string>("Thermodynamic Database Format", "simple");

  // create the appropriate chemistry object
  if (database_format == "simple") {
    chem_ = new SimpleThermoDatabase();
  } else if (database_format == "PFloTran Parsed Database") {
    // hammond stuff here....
  } else {
    // invalid database format, helpful error message and throw an error.
  }

  set_verbosity(static_cast<Verbosity>(parameter_list_.get<int>("Verbosity", 0)));

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

  // TODO: WTF? using <unsigned int> in the parameter list doesn't work...?
  beaker_parameters_.max_iterations = 
      static_cast<unsigned int>(parameter_list_.get<int>("Maximum Newton Iterations", 200));

  // physical parameters
  current_porosity_ = chemistry_state_->get_porosity();
  beaker_parameters_.porosity = (*current_porosity_)[0];

  current_water_saturation_ = chemistry_state_->get_water_saturation();
  beaker_parameters_.saturation = (*current_water_saturation_)[0];

  current_water_density_ = chemistry_state_->get_water_density();
  beaker_parameters_.water_density = (*current_water_density_)[0];

  // copy the initial conditions from some file....

  beaker_components_.free_ion.clear();
  beaker_components_.minerals.clear();
  beaker_components_.ion_exchange_sites.clear();
  beaker_components_.total.clear();
  beaker_components_.total_sorbed.clear();

  Teuchos::RCP<const Epetra_MultiVector> aqueous_components = 
      chemistry_state_->get_total_component_concentration();
  num_aqueous_components(aqueous_components->NumVectors());


  for (unsigned int c = 0; c < num_aqueous_components(); c++) {
    beaker_components_.total.push_back( 0.1 );
  }

  chem_->verbosity(verbosity());
  chem_->Setup(beaker_components_, beaker_parameters_);
  if (verbosity() > kTerse) {
    chem_->Display();
  }
  // solve for free-ion concentrations
  chem_->Speciate(beaker_components_, beaker_parameters_);
  if (verbosity() > kTerse) {
    chem_->DisplayResults();
  }
}  // end Chemistry_PK()

Chemistry_PK::~Chemistry_PK()
{
  delete chem_;
}  // end ~Chemistry_PK()


// the MPC will call this function to advance the state
// with this particular process kernel

// this is how to get the total component concentration
// CS->get_total_component_concentration()

// please update the argument to this function called tcc_star
// with the result of your chemistry computation which is
// the total component concentration ^star

// see the Chemistry_State for the other available
// data in the chemistry specific state

void Chemistry_PK::advance(
    const double& delta_time,
    Teuchos::RCP<Epetra_MultiVector> total_component_concentration_star)
{
  if (chem_->verbosity() == kDebugChemistryProcessKernel) {
    cout << "  Chemistry_PK::advance() : advancing the chemistry process model..." << endl;
  }

  current_time_ = saved_time_ + delta_time;

  // loop over all cells
  // use size of the porosity vector as indicator of size for now....
  int num_cells = chemistry_state_->get_porosity()->MyLength();
  for (int cell = 0; cell < num_cells; cell++) {

    // copy the components for this cell out of tcc or tcc_star
    for (unsigned int component = 0; component < num_aqueous_components();
         component++) {

    }

    // copy the minerals, surface complexation sites, ion exchange sites, etc



    // chemistry computations for this cell
    chem_->ReactionStep(&beaker_components_, beaker_parameters_, delta_time);

    // extract the results and save them into...? do we really want to
    // over write tcc*?
    for (unsigned int component = 0; component < num_aqueous_components();
         component++) {

    }

    // update the minerals, surface complexation sites, ion exchange sites, etc



  }
  if (chem_->verbosity() == kDebugChemistryProcessKernel) {
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
