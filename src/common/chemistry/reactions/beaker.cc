/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Driver class for evaluating geochemical related processes at a
  single computational node

  TODO(bandre): update mineral volume fractions to state.minerals after kinetics...
  TODO(bandre): finish implementing ion exchange jacobian
*/

#include <cstdlib>
#include <cassert>

#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>

#include "exceptions.hh"

#include "activity_model.hh"
#include "activity_model_factory.hh"
#include "aqueous_equilibrium_complex.hh"
#include "general_rxn.hh"
#include "radioactive_decay.hh"
#include "ion_exchange_rxn.hh"
#include "sorption_isotherm.hh"
#include "kinetic_rate.hh"
#include "mineral.hh"
#include "mineral_kinetics_factory.hh"
#include "species.hh"
#include "surface_complexation_rxn.hh"
#include "lu_solver.hh"
#include "matrix_block.hh"
#include "chemistry_utilities.hh"
#include "chemistry_exception.hh"

#include "beaker.hh"

namespace Amanzi {
namespace AmanziChemistry {

// solver defaults
const double Beaker::tolerance_default_ = 1.0e-12;
const unsigned int Beaker::max_iterations_default_ = 250;
// default physical parameters
const double Beaker::porosity_default_ = 1.0;  // [-]
const double Beaker::saturation_default_ = 1.0;  // [-]
const double Beaker::water_density_kg_m3_default_ = 1000.0;
const double Beaker::volume_default_ = 1.0;  // [m^3]

Beaker::Beaker(const Teuchos::Ptr<VerboseObject> vo)
    : vo_(vo),
      tolerance_(tolerance_default_),
      max_iterations_(max_iterations_default_),
      ncomp_(0),
      total_(),
      dtotal_(),
      total_sorbed_(),
      dtotal_sorbed_(),
      porosity_(porosity_default_),
      saturation_(saturation_default_),
      water_density_kg_m3_(water_density_kg_m3_default_),
      water_density_kg_L_(1.0),
      volume_(volume_default_),
      dt_(1.0),
      aqueous_accumulation_coef_(0.0),
      sorbed_accumulation_coef_(0.0),
      por_sat_den_vol_(0.0),
      activity_model_(NULL),
      fixed_accumulation_(),
      residual_(),
      prev_molal_(),
      rhs_(),
      jacobian_(),
      lu_solver_(),
      use_log_formulation_(true),
      sorption_isotherm_params_(4, 0.0) {
  primary_species_.clear();
  minerals_.clear();
  aqComplexRxns_.clear();
  generalKineticRxns_.clear();
  radioactive_decay_rxns_.clear();
  mineral_rates_.clear();
  ion_exchange_rxns_.clear();
  sorption_isotherm_rxns_.clear();
  total_.clear();
  total_sorbed_.clear();
}


Beaker::~Beaker() {
  delete activity_model_;

  if (mineral_rates_.size() != 0) {
    for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
         rate != mineral_rates_.end(); rate++) {
      delete(*rate);
    }
  }
}


/********************************************************************
 **
 ** public interface for drivers (batch_chem, unstructured and
 ** structured process kernels)
 **
 *******************************************************************/


/* ******************************************************************
* public setup related
****************************************************************** */
void Beaker::Setup(const BeakerState& state,
                   const BeakerParameters& parameters) {
  SetParameters(parameters);

  SetupActivityModel(parameters.activity_model_name,
                     parameters.pitzer_database, parameters.jfunction_pitzer);
  ResizeInternalMemory(static_cast<int>(primary_species().size()));
  VerifyState(state);
}


Beaker::BeakerParameters Beaker::GetDefaultParameters() const {
  Beaker::BeakerParameters parameters;

  parameters.thermo_database_file.clear();

  parameters.tolerance = tolerance_default_;
  parameters.max_iterations = max_iterations_default_;

  parameters.activity_model_name = ActivityModelFactory::unit;
  parameters.pitzer_database.clear();

  parameters.porosity = porosity_default_;
  parameters.saturation = saturation_default_;
  parameters.water_density = water_density_kg_m3_default_;  // kg / m^3
  parameters.volume = volume_default_;  // m^3

  return parameters;
}


/* ******************************************************************
* Take a parameters object that was created by the driver, and map
* the data into the appropriate chemistry object, potentially over
* riding some of our internal database data.
****************************************************************** */
void Beaker::SetParameters(const Beaker::BeakerParameters& parameters) {
  set_tolerance(parameters.tolerance);
  set_max_iterations(parameters.max_iterations);
  set_porosity(parameters.porosity);
  water_density_kg_m3(parameters.water_density);  // den = [kg/m^3]
  set_saturation(parameters.saturation);
  set_volume(parameters.volume);  // vol = [m^3]

  // calculates the coefficient in aqueous portion of accumulation term
  set_aqueous_accumulation_coef(porosity_ * saturation_ * volume_ * 1000.0 / dt_);
  set_sorbed_accumulation_coef(volume_ / dt_);

  // calculates product of porosity,saturation,water_density[kg/m^3],volume
  por_sat_den_vol(porosity_ * saturation_ * water_density_kg_m3() * volume_);
}


//
// public "computation engine" routines
//

/* ******************************************************************
* if no water density provided, default is 1000.0 kg/m^3
****************************************************************** */
int Beaker::Speciate(BeakerState* state,
                     const BeakerParameters& parameters) {
  std::stringstream message;
  double speciation_tolerance = 1.e-12;
  double residual_tolerance = 1.e-12;
  ResetStatus();
  UpdateParameters(parameters, 1.0);  // NOTE: need dt=1 to avoid divide by zero
  CheckChargeBalance(state->total);

  CopyStateToBeaker(*state);

  // store current molalities
  for (int i = 0; i < ncomp_; i++) {
    prev_molal_.at(i) = primary_species().at(i).molality();
  }

  double max_rel_change;
  int max_rel_index;
  double max_residual;
  unsigned int num_iterations = 0;
  bool calculate_activity_coefs = false;

  do {
    UpdateActivityCoefficients();
    UpdateEquilibriumChemistry();
    CalculateDTotal();

    // calculate residual
    // units of residual: mol/sec
    for (int i = 0; i < ncomp_; i++) {
      residual_.at(i) = total_.at(i) - state->total.at(i);
    }
    // add derivatives of total with respect to free to Jacobian
    // units of Jacobian: kg water/sec
    jacobian_.Zero();
    CalculateDTotal();
    jacobian_.AddValues(&dtotal_);

    rhs_ = residual_;

    // scale the Jacobian
    ScaleRHSAndJacobian();

    // for derivatives with respect to ln concentration, scale columns
    // by primary species concentrations
    for (int i = 0; i < ncomp_; i++) {
      jacobian_.ScaleColumn(i, primary_species().at(i).molality());
    }

    // call solver
    lu_solver_.Solve(&jacobian_, &rhs_);

    // calculate update truncating at a maximum of 5 in log space
    UpdateMolalitiesWithTruncation(5.0);
    // calculate maximum relative change in concentration over all species
    CalculateMaxRelChangeInMolality(&max_rel_change, &max_rel_index);

    num_iterations++;

    max_residual = 0.0;
    for (int i = 0; i < ncomp_; i++) {
      max_residual = std::max(max_residual, std::fabs(residual_.at(i)));
    }

    // if max_rel_change small enough, turn on activity coefficients
    if (max_rel_change < speciation_tolerance ||
        max_residual < residual_tolerance) {
      calculate_activity_coefs = true;
    }

    // exist if maximum relative change is below tolerance
  } while (max_rel_change > speciation_tolerance &&
           num_iterations < max_iterations_);
  
  // for now, initialize total sorbed concentrations based on the current free
  // ion concentrations
  UpdateEquilibriumChemistry();
  CopyBeakerToState(state);
  status_.num_newton_iterations = num_iterations;
  if (max_rel_change < tolerance_) {
    status_.converged = true;
  }

  return num_iterations;
}


/* ******************************************************************
* Note: the input parameter state is modified by this function.
* Initially it contains the initial component concentrations.
* On return it contains the modified values of the state.
****************************************************************** */
int Beaker::ReactionStep(BeakerState* state,
                         const BeakerParameters& parameters,
                         double dt) {
  // update class paramters
  // water_density [kg/m^3]
  // volume [m^3]
  ResetStatus();
  UpdateParameters(parameters, dt);
  CheckChargeBalance(state->total);
  CopyStateToBeaker(*state);

  // store current molalities
  for (int i = 0; i < ncomp_; i++) {
    prev_molal_.at(i) = primary_species().at(i).molality();
  }

  // initialize to a large number (not necessary, but safe)
  double max_rel_change = 1.e20;
  int max_rel_index = -1;
  unsigned int num_iterations = 0;

  // set_use_log_formulation(false);

  // lagging activity coefficients by a time step in this case
  // UpdateActivityCoefficients();

  // calculate portion of residual at time level t
  CalculateFixedAccumulation(state->total, state->total_sorbed,
                             &fixed_accumulation_);

  do {
    // update equilibrium and kinetic chemistry (rates, ion activity
    // products, etc.)
    UpdateEquilibriumChemistry();
    UpdateKineticChemistry();

    // units of residual: mol/sec
    CalculateResidual();
    // units of Jacobian: kg water/sec
    CalculateJacobian();
    // therefore, units of solution: mol/kg water (change in molality)

    rhs_ = residual_;

    // scale the Jacobian
    ScaleRHSAndJacobian();

    if (use_log_formulation_) {
      // for derivatives with respect to ln concentration, scale columns
      // by primary species concentrations  
      for (int i = 0; i < ncomp_; i++) {
        jacobian_.ScaleColumn(i, primary_species().at(i).molality());
      }
    }

    // solve J dlnc = r
    lu_solver_.Solve(&jacobian_, &rhs_);

    // units of solution: mol/kg water (change in molality)
    // calculate update truncating at a maximum of 5 in nat log space
    // update with exp(-dlnc)
    UpdateMolalitiesWithTruncation(5.0);
    // calculate maximum relative change in concentration over all species
    CalculateMaxRelChangeInMolality(&max_rel_change, &max_rel_index);

    num_iterations++;

    //  if (num_iterations >= 100) {
    //    for (int i = 0; i < ncomp_; i++)
    //      std::cout << primary_species_.at(i).name() << " " <<
    //                   primary_species_.at(i).molality() << " " << total_.at(i) << "\n";
    //      std::cout << max_rel_change << " " << tolerance_ << std::endl;
    //  }

    // exit if maximum relative change is below tolerance
  } while (max_rel_change > tolerance_ && num_iterations < max_iterations_);

  if (num_iterations >= max_iterations_) {
    // TODO(bandre): should this be an error to the driver...?
    // code eventually produces nans when this isn't an error.
    std::ostringstream error_stream;
    error_stream << "Warning: The maximum number Netwon iterations reached in Beaker::ReactionStep()." << std::endl;
    error_stream << "Warning: Results may not have the desired accuracy." << std::endl;
    error_stream << "Warning: max relative change = " << max_rel_change << std::endl;
    error_stream << "Warning: max relative index = " << max_rel_index << std::endl;
    error_stream << "Warning: tolerance = " << tolerance_ << std::endl;
    error_stream << "Warning: max iterations = " << max_iterations_ << std::endl;
    // update before leaving so that we can see the erroneous values!
    CopyBeakerToState(state);
    Exceptions::amanzi_throw(ChemistryMaxIterationsReached(error_stream.str()));
  }

  status_.num_newton_iterations = num_iterations;
  if (max_rel_change < tolerance_) {
    status_.converged = true;
  }

  // update total concentrations
  UpdateEquilibriumChemistry();
  UpdateKineticMinerals();

  // TODO(bandre): not convinced this is the correct place to call
  // UpdateActivityCoefficients. Should be before the call to
  // UpdateEquilibriumChemistry()? But that changes the numerical
  // results and I need to look more closely at what is going on.

  // lagging activity coefficients by a time step
  UpdateActivityCoefficients();

  CopyBeakerToState(state);
  ValidateSolution();

  return num_iterations;
}



/* ******************************************************************
* Copy the beaker state into variables are be returned to the
* driver for long term storage.
****************************************************************** */
void Beaker::CopyBeakerToState(Beaker::BeakerState* state) {
  // NOTE: The state struct may have been only partially
  // initialized be the driver (it doesn't know about the size of
  // internal memory for surface complexation or ion exchange....
  // For now we assert the things the driver should know....

  // totals
  assert(state->total.size() == total_.size());
  for (int i = 0; i < ncomp_; ++i) {
    state->total.at(i) = total_.at(i);
  }

  if (total_sorbed_.size() > 0) {
    assert(state->total_sorbed.size() == total_sorbed_.size());
    for (int i = 0; i < ncomp_; ++i) {
      state->total_sorbed.at(i) = total_sorbed_.at(i);
    }
  }

  // free ion
  assert(state->free_ion.size() == ncomp_);
  for (int i = 0; i < ncomp_; ++i) {
    state->free_ion.at(i) = primary_species().at(i).molality();
  }

  //
  // activity coeff
  //
  if (state->primary_activity_coeff.size() != ncomp_) {
    state->primary_activity_coeff.resize(ncomp_);
  }
  for (int i = 0; i < ncomp_; ++i) {
    state->primary_activity_coeff.at(i) = primary_species().at(i).act_coef();
  }
  if (state->secondary_activity_coeff.size() != aqComplexRxns_.size()) {
    state->secondary_activity_coeff.resize(aqComplexRxns_.size());
  }
  for (int i = 0; i < aqComplexRxns_.size(); ++i) {
    state->secondary_activity_coeff.at(i) = aqComplexRxns_.at(i).act_coef();
  }

  //
  // minerals
  //
  assert(state->mineral_volume_fraction.size() == minerals_.size());
  for (unsigned int m = 0; m < minerals_.size(); ++m) {
    state->mineral_volume_fraction.at(m) = minerals_.at(m).volume_fraction();
  }

  if (state->mineral_specific_surface_area.size() != minerals_.size()) {
    state->mineral_specific_surface_area.resize(minerals_.size());
  }
  for (int m = 0; m < minerals_.size(); ++m) {
    double ssa = minerals_.at(m).specific_surface_area();
    state->mineral_specific_surface_area.at(m) = ssa;
  }

  //
  // ion exchange
  //
  if (state->ion_exchange_sites.size() != ion_exchange_rxns_.size()) {
    state->ion_exchange_sites.resize(ion_exchange_rxns_.size());
  }
  for (int i = 0; i < ion_exchange_rxns_.size(); ++i) {
    state->ion_exchange_sites.at(i) = 
        ion_exchange_rxns_.at(i).site().get_cation_exchange_capacity();
  }

  if (state->ion_exchange_ref_cation_conc.size() != ion_exchange_rxns_.size()) {
    state->ion_exchange_ref_cation_conc.resize( 
        ion_exchange_rxns_.size());
  }
  for (int i = 0; i < ion_exchange_rxns_.size(); ++i) {
    state->ion_exchange_ref_cation_conc.at(i) = 
        ion_exchange_rxns_[i].ref_cation_sorbed_conc();
  }

  //
  // surface complexation
  //
  if (state->surface_complex_free_site_conc.size() != 
      surfaceComplexationRxns_.size()) {
    state->surface_complex_free_site_conc.resize( 
        surfaceComplexationRxns_.size());
  }
  for (unsigned int i = 0; i < surfaceComplexationRxns_.size(); ++i) {
    state->surface_complex_free_site_conc.at(i) =
        surfaceComplexationRxns_.at(i).free_site_concentration();
  }

  if (state->surface_site_density.size() != 
      surfaceComplexationRxns_.size()) {
    state->surface_site_density.resize(surfaceComplexationRxns_.size());
  }
  for (unsigned int i = 0; i < surfaceComplexationRxns_.size(); ++i) {
    state->surface_site_density.at(i) =
        surfaceComplexationRxns_.at(i).GetSiteDensity();
  }

  //
  // sorption isotherms
  //
  if (sorption_isotherm_rxns_.size() > 0) {
    if (state->isotherm_kd.size() != ncomp_) {
      state->isotherm_kd.resize(ncomp_, 0.0);
    }
    if (state->isotherm_langmuir_b.size() != ncomp_) {
      state->isotherm_langmuir_b.resize(ncomp_, 0.0);
    }
    if (state->isotherm_freundlich_n.size() != ncomp_) {
      state->isotherm_freundlich_n.resize(ncomp_, 1.0);
    }
    for (int r = 0; r < sorption_isotherm_rxns_.size(); ++r) {
      const std::vector<double>& params = sorption_isotherm_rxns_.at(r).GetIsothermParameters();
      int id = sorption_isotherm_rxns_.at(r).species_id();
      state->isotherm_kd.at(id) = params.at(0);
      if (sorption_isotherm_rxns_.at(r).IsothermType() == SorptionIsotherm::FREUNDLICH) {
        state->isotherm_freundlich_n.at(id) = params.at(1);
      } else if (sorption_isotherm_rxns_.at(r).IsothermType() == SorptionIsotherm::LANGMUIR) {
        state->isotherm_langmuir_b.at(id) = params.at(1);
      }
    }
  }
}


void Beaker::GetPrimaryNames(std::vector<std::string>* names) const {
  names->clear();
  for (std::vector<Species>::const_iterator primary = primary_species().begin();
       primary != primary_species().end(); primary++) {
    names->push_back(primary->name());
  }
}


int Beaker::GetPrimaryIndex(const std::string& name) const {
  int index = -1;
  for (std::vector<Species>::const_iterator primary = primary_species().begin();
       primary != primary_species().end(); primary++) {
    if (primary->name() == name) {
      index = primary->identifier();
    }
  }
  return index;
}


//
// public output
//

void Beaker::Display() const {
  vo_->Write(Teuchos::VERB_HIGH, "-- Beaker description ------------------------------------------------\n");
  DisplayParameters();

  DisplayPrimary();

  DisplayAqueousEquilibriumComplexes();

  DisplayGeneralKinetics();

  DisplayRadioactiveDecayRxns();

  DisplayMinerals();

  DisplayMineralKinetics();

  DisplayIonExchangeSites();

  DisplayIonExchangeComplexes();

  DisplaySurfaceSites();

  DisplaySurfaceComplexes();

  DisplaySorptionIsotherms();

  vo_->Write(Teuchos::VERB_HIGH, "------------------------------------------------ Beaker description --\n");
}


void Beaker::DisplayComponents(const Beaker::BeakerState& state) const {
  std::stringstream message;
  message << "--- Input Components -------------------------------------------------"
          << std::endl;
  message << "---- Aqueous Components" << std::endl;
  message << std::setw(15) << "Name"
          << std::setw(15) << "Molality"
          << std::setw(15) << "Molarity"
          // << std::setw(15) << "Free Ion" // TODO(bandre): uncomment and update test results
          << std::endl;
  for (int i = 0; i < ncomp_; i++) {
    message << std::setw(15) << primary_species().at(i).name()
            << std::scientific << std::setprecision(5)
            << std::setw(15) << state.total.at(i) / water_density_kg_L()
            << std::setw(15) << state.total.at(i)
            // << std::setw(15) << state.free_ion.at(i) // TODO(bandre): uncomment and update test results
            << std::endl;
  }

  if (minerals_.size() > 0) {
    message << "---- Mineral Components" << std::endl;
    message << std::setw(15) << "Name"
            << std::setw(15) << "Vol. frac" << std::endl;
    for (unsigned int m = 0; m < minerals_.size(); m++) {
      message << std::setw(15) << minerals_.at(m).name()
              << std::setw(15) << std::fixed << std::setprecision(5)
              << state.mineral_volume_fraction.at(m) << std::endl;
    }
  }

  if (total_sorbed_.size() > 0) {
    message << "---- Sorbed Components" << std::endl;
    message << std::setw(15) << "Name"
            << std::setw(15) << "Moles / m^3" << std::endl;
    for (int i = 0; i < ncomp_; i++) {
      message << std::setw(15) << primary_species().at(i).name()
              << std::scientific << std::setprecision(5)
              << std::setw(15) << state.total_sorbed.at(i)
              << std::endl;
    }
  }
  message << "------------------------------------------------- Input Components ---"
          << std::endl;
  vo_->Write(Teuchos::VERB_HIGH, message.str());
}


void Beaker::DisplayResults() const {
  std::stringstream message;
  message << std::endl;
  message << "-- Solution ----------------------------------------------------------"
          << std::endl;
  message << "---- Components " << std::endl;
  message << std::setw(15) << "Name"
          << std::setw(15) << "Molality"
          << std::setw(15) << "Molarity"
          << std::endl;
  for (int i = 0; i < ncomp_; i++) {
    message << std::setw(15) << primary_species().at(i).name()
            << std::scientific << std::setprecision(5)
            << std::setw(15) << total_.at(i) / water_density_kg_L()
            << std::setw(15) << total_.at(i)
            << std::endl;
  }

  message << "---- Change Balance " << std::endl;
  double charge_balance_molal = 0.0;
  for (int i = 0; i < ncomp_; i++) {
    charge_balance_molal += primary_species().at(i).charge() * total_.at(i);
  }
  message << std::setw(15) << " "
          << std::scientific << std::setprecision(5)
          << std::setw(15) << " "
          << std::setw(15) << charge_balance_molal
          << std::endl;

  message << "---- Species " << std::endl;
  vo_->Write(Teuchos::VERB_HIGH, message.str());

  primary_species().at(0).DisplayResultsHeader(vo_);
  for (int i = 0; i < ncomp_; i++) {
    primary_species().at(i).DisplayResults(vo_);
  }

  // same header info as primaries....
  for (unsigned int i = 0; i < aqComplexRxns_.size(); i++) {
    aqComplexRxns_.at(i).DisplayResults(vo_);
  }

  if (minerals_.size() > 0) {
    vo_->Write(Teuchos::VERB_HIGH, "---- Minerals\n");
    minerals_[0].DisplayResultsHeader(vo_);
    for (unsigned int i = 0; i < minerals_.size(); i++) {
      minerals_.at(i).DisplayResults(vo_);
    }
  }

  if (ion_exchange_rxns_.size() > 0) {
    vo_->Write(Teuchos::VERB_HIGH, "---- Ion Exchange Sites\n");
    ion_exchange_rxns_.at(0).site().DisplayResultsHeader(vo_);
    for (unsigned int i = 0; i < ion_exchange_rxns_.size(); i++) {
      ion_exchange_rxns_.at(i).site().DisplayResults(vo_);
    }
  }

  if (ion_exchange_rxns_.size() > 0) {
    vo_->Write(Teuchos::VERB_HIGH, "---- Ion Exchange Complexes\n");
    ion_exchange_rxns_.at(0).ionx_complexes().at(0).DisplayResultsHeader(vo_);
    for (unsigned int i = 0; i < ion_exchange_rxns_.size(); i++) {
      for (unsigned int j = 0; j < ion_exchange_rxns_.at(i).ionx_complexes().size(); j++) {
        (ion_exchange_rxns_.at(i).ionx_complexes())[j].DisplayResults(vo_);
      }
    }
  }

  if (surfaceComplexationRxns_.size() > 0) {
    vo_->Write(Teuchos::VERB_HIGH, "---- Surface Complexation Reactions\n");
    for (unsigned int i = 0; i < surfaceComplexationRxns_.size(); i++) {
      surfaceComplexationRxns_.at(i).DisplayResultsHeader(vo_);
      surfaceComplexationRxns_.at(i).DisplayResults(vo_);
    }
  }

  vo_->Write(Teuchos::VERB_HIGH, "---------------------------------------------------------- Solution --\n\n");
}


void Beaker::DisplayTotalColumnHeaders(const bool display_free) const {
  std::stringstream message;
  message << std::setw(15) << "Time (s)";
  for (int i = 0; i < ncomp_; i++) {
    message << std::setw(15) << primary_species().at(i).name();
  }
  if (display_free) {
    for (int i = 0; i < ncomp_; i++) {
      std::string temp = primary_species().at(i).name() + "_free";
      message << std::setw(15) << temp;
    }
  }
  if (total_sorbed_.size() > 0) {
    for (int i = 0; i < total_sorbed_.size(); i++) {
      std::string temp = primary_species().at(i).name() + "_sorbed";
      message << std::setw(15) << temp;
    }
  }
  if (minerals().size() > 0) {
    for (int m = 0; m < minerals().size(); ++m) {
      std::string temp =  minerals().at(m).name() + "_vf";
      message << std::setw(15) << temp;
    }
  }
  message << std::endl;
  vo_->Write(Teuchos::VERB_HIGH, message.str());
}


void Beaker::DisplayTotalColumns(const double time,
                                 const BeakerState& state,
                                 const bool display_free) const {
  std::stringstream message;
  message << std::scientific << std::setprecision(6) << std::setw(15);
  message << time;
  for (int i = 0; i < ncomp_; i++) {
    message << std::setw(15) << state.total.at(i);
  }
  if (display_free) {
    for (int i = 0; i < ncomp_; i++) {
      message << std::setw(15) << state.free_ion.at(i);
    }
  }
  if (total_sorbed_.size() > 0) {
    for (int i = 0; i < total_sorbed_.size(); i++) {
      message << std::setw(15) << state.total_sorbed.at(i);
    }
  }
  if (minerals().size() > 0) {
    for (int m = 0; m < minerals().size(); ++m) {
      message << std::setw(15) << state.mineral_volume_fraction.at(m);
    }
  }
  message << std::endl;
  vo_->Write(Teuchos::VERB_HIGH, message.str());
}


/*******************************************************************************
 **
 ** protected interface (SimpleThermoDatabase)
 **
 ******************************************************************************/
void Beaker::ResizeInternalMemory(const int size) {
  set_ncomp(size);
  total_.resize(size);
  dtotal_.Resize(size);

  if (surfaceComplexationRxns_.size() > 0 ||
      sorption_isotherm_rxns_.size() > 0 ||
      ion_exchange_rxns_.size() > 0) {
    total_sorbed_.resize(size, 0.0);
    dtotal_sorbed_.Resize(size);
    dtotal_sorbed_.Zero();
  } else {
    total_sorbed_.resize(0);
    // dtotal_sorbed_.Resize(0);
  }

  fixed_accumulation_.resize(size);
  residual_.resize(size);
  prev_molal_.resize(size);

  jacobian_.Resize(size);
  rhs_.resize(size);
  lu_solver_.Initialize(size);
}


/* ******************************************************************
* some helpful error checking goes here...
****************************************************************** */
void Beaker::VerifyState(const Beaker::BeakerState& state) const {
  bool error = false;
  std::ostringstream error_stream;
  error_stream << "Beaker::VerifyState():\ndatabase input and component initial conditions do not match:\n";

  // if the size of the various initial conditions, state, and
  // database input don't match. Print a helpful message and exit
  // gracefully.

  if (static_cast<unsigned int>(ncomp_) != state.total.size()) {
    error = true;
    error_stream << "ncomp(" << ncomp_
                 << ") and state.total.size(" << state.total.size()
                 << ") do not match.\n";
  }

  if (primary_species().size() != state.total.size()) {
    error = true;
    error_stream << "primary_species.size(" << this->primary_species().size()
                 << ") and state.total.size(" << state.total.size()
                 << ") do not match.\n";
  }

  if (state.free_ion.size() != state.total.size()) {
    error = true;
    error_stream << "state.total.size(" << state.total.size()
                 << ") and state.free_ion.size(" << state.free_ion.size()
                 << ") do not match.\n";
  }

  if (ion_exchange_rxns().size() != state.ion_exchange_sites.size()) {
    error = true;
    error_stream << "ion_exchange_rxns.size(" << this->ion_exchange_rxns().size()
                 << ") and state.ion_exchange_sites.size(" << state.ion_exchange_sites.size()
                 << ") do not match.\n";
  }

  if (minerals().size() != state.mineral_volume_fraction.size()) {
    error = true;
    error_stream << "minerals.size(" << this->minerals().size()
                 << ") and state.mineral_volume_fraction.size(" << state.mineral_volume_fraction.size()
                 << ") do not match.\n";
  }

  if (total().size() != state.total.size()) {
    error = true;
    error_stream << "total.size(" << this->total().size()
                 << ") and state.total.size(" << state.total.size()
                 << ") do not match.\n";
  }

  // FIXED? this check is breaking things because total_sorbed is always
  // resized in resize(), even if there is no sorption!
  if (total_sorbed().size() != state.total_sorbed.size()) {
    error = true;
    error_stream << "total_sorbed.size(" << this->total_sorbed().size()
                 << ") and state.total_sorbed.size("
                 << state.total_sorbed.size() << ") do not match.\n";
  }

  if (error) {
    Exceptions::amanzi_throw(ChemistryMemorySizeError(error_stream.str()));
  }
}


/* ******************************************************************
* NOTE: Do not copy total and total_sorbed here!
****************************************************************** */
void Beaker::CopyStateToBeaker(const Beaker::BeakerState& state) {
  // free ion
  if (state.free_ion.size() > 0) {
    InitializeMolalities(state.free_ion);
  } else {
    InitializeMolalities(1.e-9);
  }

  // activity coefficients
  if (state.primary_activity_coeff.size() > 0) {
    assert(state.primary_activity_coeff.size() == primary_species().size());
    for (int i = 0; i < primary_species().size(); ++i) {
      double value = state.primary_activity_coeff.at(i);
      primary_species_.at(i).set_act_coef(value);
    }
  }
  if (state.secondary_activity_coeff.size() > 0) {
    assert(state.secondary_activity_coeff.size() == aqComplexRxns_.size());
    for (int i = 0; i < aqComplexRxns_.size(); ++i) {
      double value = state.secondary_activity_coeff.at(i);
      aqComplexRxns_.at(i).set_act_coef(value);
    }
  }

  // minerals
  unsigned int size = state.mineral_volume_fraction.size();
  if (minerals().size() == size) {
    for (unsigned int m = 0; m < size; m++) {
      minerals_.at(m).set_volume_fraction(state.mineral_volume_fraction.at(m));
    }
  } else {
    std::ostringstream error_stream;
    error_stream << "Beaker::CopyStateToBeaker(): \n";
    error_stream << "minerals.size(" << minerals().size()
                 << ") and state.mineral_volume_fraction.size(" << size
                 << ") do not match.\n";
    Exceptions::amanzi_throw(ChemistryUnrecoverableError(error_stream.str()));
  }

  if (state.mineral_specific_surface_area.size() > 0) {
    assert(state.mineral_specific_surface_area.size() == minerals_.size());
    for (unsigned int m = 0; m < minerals_.size(); ++m) {
      double value = state.mineral_specific_surface_area.at(m);
      minerals_.at(m).set_specific_surface_area(value);
    }
  }

  // ion exchange
  size = state.ion_exchange_sites.size();
  if (ion_exchange_rxns().size() == size) {
    for (unsigned int ies = 0; ies < size; ies++) {
      ion_exchange_rxns_[ies].set_cation_exchange_capacity(state.ion_exchange_sites.at(ies));
    }
  } else {
    std::ostringstream error_stream;
    error_stream << "Beaker::CopyStateToBeaker(): \n";
    error_stream << "ion_exchange_sites.size(" << ion_exchange_rxns().size()
                 << ") and state.ion_exchange_sites.size(" << size 
                 << ") do not match.\n";
    Exceptions::amanzi_throw(ChemistryUnrecoverableError(error_stream.str()));
  }

  if (ion_exchange_rxns().size() > 0) {
    // check for the ref cation conc....
    if (state.ion_exchange_ref_cation_conc.size() ==
        ion_exchange_rxns().size()) {
      // we have a value from a previous solve, restore it.
      for (int i = 0; i < ion_exchange_rxns().size(); ++i) {
        double value = state.ion_exchange_ref_cation_conc.at(i);
        ion_exchange_rxns_[i].set_ref_cation_sorbed_conc(value);
      }
    } else {
      // no previous value, provide a guess
      for (int i = 0; i < ion_exchange_rxns().size(); ++i) {
        ion_exchange_rxns_[i].set_ref_cation_sorbed_conc(1.0e-9);
      }
    }
  }

  // surface complexation
  if (state.surface_site_density.size() > 0) {
    assert(state.surface_site_density.size() == surfaceComplexationRxns_.size());
    for (unsigned int i = 0; i < surfaceComplexationRxns_.size(); ++i) {
      double value = state.surface_site_density.at(i);
      surfaceComplexationRxns_.at(i).UpdateSiteDensity(value);
    }
  }

  if (surfaceComplexationRxns_.size() > 0) {
    if (state.surface_complex_free_site_conc.size() ==
        surfaceComplexationRxns_.size()) {
      // we have a value from a previous solve, restore it.
      for (int r = 0; r < surfaceComplexationRxns_.size(); ++r) {
        double value = state.surface_complex_free_site_conc.at(r);
        surfaceComplexationRxns_.at(r).set_free_site_concentration(value);
      }
    } else {
      // no previous value, provide a guess
      for (int r = 0; r < surfaceComplexationRxns_.size(); ++r) {
        //double value = 0.1 * surfaceComplexationRxns_.at(r).GetSiteDensity();
        double value = 1.0e-9;
        surfaceComplexationRxns_.at(r).set_free_site_concentration(value);
      }
    }
  }

  // sorption isotherms
  if (sorption_isotherm_rxns_.size() > 0 &&
      state.isotherm_kd.size() > 0) {
    // the driver maybe attempting to over the database values, or the
    // state were resized by a call to CopyBeakerToState()
    assert(state.isotherm_kd.size() == ncomp_);
    assert(state.isotherm_freundlich_n.size() == ncomp_);
    assert(state.isotherm_langmuir_b.size() == ncomp_);
    // NOTE(bandre): sorption_isotherm_params_ hard coded size=4,
    // current max parameters is 2
    for (int r = 0; r < sorption_isotherm_rxns_.size(); ++r) {
      int id = sorption_isotherm_rxns_.at(r).species_id();
      sorption_isotherm_params_.at(0) = state.isotherm_kd.at(id);
      SorptionIsotherm::SorptionIsothermType isotherm_type = sorption_isotherm_rxns_.at(r).IsothermType();
      if (isotherm_type == SorptionIsotherm::FREUNDLICH) {
        sorption_isotherm_params_.at(1) = state.isotherm_freundlich_n.at(id);
      } else if (isotherm_type == SorptionIsotherm::LANGMUIR) {
        sorption_isotherm_params_.at(1) = state.isotherm_langmuir_b.at(id);
      }
      sorption_isotherm_rxns_.at(r).SetIsothermParameters(sorption_isotherm_params_);
    }
  }
} 


void Beaker::SetupActivityModel(std::string model,
                                std::string pitzer_database,
                                std::string jfunction_pitzer) {

  delete activity_model_;

  ActivityModel::ActivityModelParameters parameters;
  parameters.database_filename = pitzer_database;
  parameters.pitzer_jfunction = jfunction_pitzer;

  ActivityModelFactory amf;

  activity_model_ = amf.Create(model, parameters,
                               primary_species(), aqComplexRxns_,
                               vo_);

  // activity_model_->Display();
}


void Beaker::AddPrimarySpecies(const Species& s) {
  primary_species_.push_back(s);
}


void Beaker::AddIonExchangeRxn(const IonExchangeRxn& ionx_rxn) {
  ion_exchange_rxns_.push_back(ionx_rxn);
}


void Beaker::AddIonExchangeComplex(const int irxn, const IonExchangeComplex& ionx_complex) {
  ion_exchange_rxns_[irxn].AddIonExchangeComplex(ionx_complex);
}


void Beaker::AddAqueousEquilibriumComplex(const AqueousEquilibriumComplex& c) {
  aqComplexRxns_.push_back(c);
}


void Beaker::AddMineral(const Mineral& m) {
  minerals_.push_back(m);
}


void Beaker::AddMineralKineticRate(KineticRate* rate) {
  mineral_rates_.push_back(rate);
}


bool Beaker::HaveKinetics() const {
  bool have_kinetics = false;
  if (mineral_rates_.size()) {
    have_kinetics = true;
  }
  // add other kinetic processes here....

  return have_kinetics;
}


void Beaker::AddGeneralRxn(const GeneralRxn& r) {
  generalKineticRxns_.push_back(r);
}


void Beaker::AddRadioactiveDecayRxn(const RadioactiveDecay& r) {
  radioactive_decay_rxns_.push_back(r);
}


void Beaker::AddSurfaceComplexationRxn(const SurfaceComplexationRxn& r) {
  surfaceComplexationRxns_.push_back(r);
}


void Beaker::AddSorptionIsothermRxn(const SorptionIsothermRxn& r) {
  sorption_isotherm_rxns_.push_back(r);
}


/*******************************************************************************
 **
 ** private routines
 **
 ******************************************************************************/
void Beaker::UpdateParameters(const Beaker::BeakerParameters& parameters,
                              double delta_t) {
  set_dt(delta_t);  // delta time = [sec]
  SetParameters(parameters);
} 


void Beaker::ResetStatus() {
  status_.num_rhs_evaluations = 0;
  status_.num_jacobian_evaluations = 0;
  status_.num_newton_iterations = 0;
  status_.converged = false;
}


void Beaker::UpdateActivityCoefficients() {
  activity_model_->CalculateIonicStrength(primary_species_, aqComplexRxns_);
  activity_model_->CalculateActivityCoefficients(&primary_species_,
                                                 &aqComplexRxns_,
                                                 &water_);
  for (auto it = primary_species_.begin(); it != primary_species_.end(); ++it) {
    it->update();
  }
}


void Beaker::UpdateKineticMinerals() {
  // TODO(bandre): Need to move this into the N-R loop and cut the
  // reaction rate or time step if volume fractions go negative. Right
  // now we are just setting volume fraction to zero and introducing
  // mass balance errors!

  // loop through the kinetic minerals list. Update the volume
  // fraction, specific surface area, etc

  for (std::vector<KineticRate*>::iterator r = mineral_rates_.begin();
       r != mineral_rates_.end(); ++r) {
    double kinetic_rate = (*r)->reaction_rate();
    int i = (*r)->identifier();
    minerals_.at(i).UpdateVolumeFraction(kinetic_rate, dt_);
    minerals_.at(i).UpdateSpecificSurfaceArea();
  }
}


void Beaker::InitializeMolalities(double initial_molality) {
  for (auto it = primary_species_.begin(); it != primary_species_.end(); ++it) {
    it->update(initial_molality);
  }
}


void Beaker::InitializeMolalities(const std::vector<double>& initial_molalities) {
  if (initial_molalities.size() != primary_species().size()) {
    std::ostringstream error_stream;
    error_stream << "Beaker::InitializeMolalities(): \n";
    error_stream << "   Mismatch in size of initial_molalities array (" 
                 << initial_molalities.size() << ") and the number of "
                 << "primary species (" << primary_species().size() << ")\n";
    Exceptions::amanzi_throw(ChemistryInvalidInput(error_stream.str()));
  }

  // iterator doesnt seem to work then passing a vector entry - geh
  for (unsigned int i = 0; i < primary_species().size(); i++) {
    primary_species_.at(i).update(initial_molalities.at(i));
  }
}


void Beaker::UpdateEquilibriumChemistry() {
  // calculateActivityCoefficients(-1);

  // update primary species activities
  for (auto it = primary_species_.begin(); it != primary_species_.end(); ++it) {
    it->update();
  }

  // calculated secondary aqueous complex concentrations
  for (auto it = aqComplexRxns_.begin(); it != aqComplexRxns_.end(); ++it) {
    it->Update(primary_species(), water_);
  }

  // calculate mineral saturation states
  for (auto it = minerals_.begin(); it != minerals_.end(); ++it) {
    it->Update(primary_species(), water_);
  }

  // surface complexation
  for (auto it = surfaceComplexationRxns_.begin(); it != surfaceComplexationRxns_.end(); ++it) {
    it->Update(primary_species());
  }

  // sorption isotherms
  for (auto it = sorption_isotherm_rxns_.begin(); it != sorption_isotherm_rxns_.end(); ++it) {
    it->Update(primary_species());
  }

  // add equilibrium ion exchange here?
  for (auto it = ion_exchange_rxns_.begin(); it != ion_exchange_rxns_.end(); ++it) {
    it->Update(primary_species());
  }

  // calculate total component concentrations
  CalculateTotal();
}


void Beaker::CalculateTotal() {
  // add in primaries
  for (unsigned int i = 0; i < total_.size(); i++) {
    total_.at(i) = primary_species().at(i).molality();
  }

  // add in aqueous complexes
  for (auto it = aqComplexRxns_.begin(); it != aqComplexRxns_.end(); ++it) {
    it->AddContributionToTotal(&total_);
  }

  // scale by water density to convert to molarity
  for (unsigned int i = 0; i < total_.size(); i++) {
    total_.at(i) *= water_density_kg_L();
  }

  // calculate sorbed totals
  // initialize to zero
  for (unsigned int i = 0; i < total_sorbed_.size(); i++) {
    total_sorbed_.at(i) = 0.0;
  }

  // add in surface complex contributions
  for (auto it = surfaceComplexationRxns_.begin(); it != surfaceComplexationRxns_.end(); ++it) {
    it->AddContributionToTotal(&total_sorbed_);
  }

  // add in isotherm contributions
  for (auto it = sorption_isotherm_rxns_.begin(); it != sorption_isotherm_rxns_.end(); ++it) {
    it->AddContributionToTotal(&total_sorbed_);
  }

  // add ion exchange
  for (auto it = ion_exchange_rxns_.begin(); it != ion_exchange_rxns_.end(); ++it) {
    it->AddContributionToTotal(&total_sorbed_);
  }
}


void Beaker::CalculateDTotal() {
  dtotal_.Zero();
  // derivative with respect to free-ion is 1.
  dtotal_.SetDiagonal(1.0);

  // add in derviative of complex contribution with respect to free-ion
  for (auto it = aqComplexRxns_.begin(); it != aqComplexRxns_.end(); ++it) {
    it->AddContributionToDTotal(primary_species(), &dtotal_);
  }

  // scale by density of water
  dtotal_.Scale(water_density_kg_L());

  // calculate sorbed derivatives
  if (total_sorbed_.size() > 0) {
    dtotal_sorbed_.Zero();
    for (auto it = surfaceComplexationRxns_.begin(); it != surfaceComplexationRxns_.end(); ++it) {
      it->AddContributionToDTotal(primary_species(), &dtotal_sorbed_);
    }
    for (auto it = sorption_isotherm_rxns_.begin(); it != sorption_isotherm_rxns_.end(); ++it) {
      it->AddContributionToDTotal(primary_species(), &dtotal_sorbed_);
    }
    // add ion exchange
    for (auto it = ion_exchange_rxns_.begin(); it != ion_exchange_rxns_.end(); ++it) {
      it->AddContributionToDTotal(primary_species(), &dtotal_sorbed_);
    }
  }
}


void Beaker::UpdateKineticChemistry() {
  // loop over general kinetic reactions and update effective rates
  for (std::vector<GeneralRxn>::iterator i = generalKineticRxns_.begin();
       i != generalKineticRxns_.end(); i++) {
    i->update_rates(primary_species());
  }

  // loop over radioactive decay reactions and update effective rates
  // NOTE(bandre): radio active decay operates on total, not free conc.
  // need to pass the volume of liquid: porosity * saturation * volume
  for (std::vector<RadioactiveDecay>::iterator i = radioactive_decay_rxns_.begin();
       i != radioactive_decay_rxns_.end(); ++i) {
    i->UpdateRate(total_, total_sorbed_, porosity_, saturation_, volume_);
  }

  // add mineral saturation and rate calculations here
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->Update(primary_species(), minerals_);
  }
  // add multirate kinetic surface complexation reaction quotient calculations
  // here
} 


void Beaker::AddKineticChemistryToResidual() {
  // loop over general kinetic reactions and add rates
  for (std::vector<GeneralRxn>::iterator i = generalKineticRxns_.begin();
       i != generalKineticRxns_.end(); i++) {
    i->addContributionToResidual(&residual_, por_sat_den_vol());
  }

  // loop over radioactive decay reactions and add rates
  for (std::vector<RadioactiveDecay>::iterator i = radioactive_decay_rxns_.begin();
       i != radioactive_decay_rxns_.end(); ++i) {
    i->AddContributionToResidual(&residual_);
  }

  // add mineral mineral contribution to residual here.  units = mol/sec.
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->AddContributionToResidual(minerals_, volume_, &residual_);
  }

  // add multirate kinetic surface complexation contribution to residual here.
}


void Beaker::AddKineticChemistryToJacobian() {
  // loop over general kinetic reactions and add rates
  for (std::vector<GeneralRxn>::iterator i = generalKineticRxns_.begin();
       i != generalKineticRxns_.end(); i++) {
    i->addContributionToJacobian(&jacobian_, primary_species(), por_sat_den_vol());
  }

  // loop over radioactive decay reactions and add rates
  for (std::vector<RadioactiveDecay>::iterator i = radioactive_decay_rxns_.begin();
       i != radioactive_decay_rxns_.end(); ++i) {
    i->AddContributionToJacobian(dtotal_, dtotal_sorbed_,
                                 porosity_, saturation_, volume_,
                                 &jacobian_);
  }

  // add mineral mineral contribution to Jacobian here.  units = kg water/sec.
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->AddContributionToJacobian(primary_species(), minerals_, volume_, &jacobian_);
  }

  // add multirate kinetic surface complexation contribution to Jacobian here.
}


void Beaker::AddAccumulation(const std::vector<double>& total,
                             const std::vector<double>& total_sorbed,
                             std::vector<double>* residual) {
  // aqueous_accumulation_coef = porosity*saturation*volume*1000./dt
  // units = (mol solute/L water)*(m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  //         (m^3 bulk)*(1000L water/m^3 water)/(sec) = mol/sec
  // 1000.d0 converts vol from m^3 -> L
  // all residual entries should be in mol/sec
  for (unsigned int i = 0; i < total.size(); i++) {
    residual->at(i) += aqueous_accumulation_coef_ * total.at(i);
  }

  // add accumulation term for equilibrium sorption (e.g. Kd, surface
  // complexation) here
  // sorbed_accumulation_coef = volume/dt
  // units = (mol solute/m^3 bulk)*(m^3 bulk)/(sec) = mol/sec
  // all residual entries should be in mol/sec
  for (unsigned int i = 0; i < total_sorbed.size(); i++) {
    residual->at(i) += sorbed_accumulation_coef() * total_sorbed.at(i);
  }
}


void Beaker::AddAccumulationDerivative(MatrixBlock* J,
                                       MatrixBlock* dtotal,
                                       MatrixBlock* dtotal_sorbed) {
  // aqueous_accumulation_coef = porosity*saturation*volume*1000./dt
  // units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*(m^3 bulk)/(sec)
  //         *(kg water/L water)*(1000L water/m^3 water) = kg water/sec
  // all Jacobian entries should be in kg water/sec
  J->AddValues(dtotal, aqueous_accumulation_coef_);

  // add accumulation derivative term for equilibrium sorption
  // (e.g. Kd, surface complexation) here
  // sorbed_accumulation_coef = volume/dt
  // units = (kg water/m^3 bulk)*(m^3 bulk)/(sec) = kg water/sec
  if (total_sorbed_.size()) {
    J->AddValues(dtotal_sorbed, sorbed_accumulation_coef());
  }
}


void Beaker::CalculateFixedAccumulation(const std::vector<double>& total,
                                        const std::vector<double>& total_sorbed,
                                        std::vector<double>* fixed_accumulation) {
  for (unsigned int i = 0; i < total.size(); i++) {
    fixed_accumulation->at(i) = 0.0;
  }
  AddAccumulation(total, total_sorbed, fixed_accumulation);
}


void Beaker::CalculateResidual() {
  status_.num_rhs_evaluations++;
  // subtract fixed porition
  for (int i = 0; i < ncomp_; i++) {
    residual_.at(i) = -fixed_accumulation_.at(i);
  }

  // accumulation adds in equilibrium chemistry
  AddAccumulation(total_, total_sorbed_, &residual_);

  // kinetic reaction contribution to residual
  AddKineticChemistryToResidual();
}


void Beaker::CalculateJacobian() {
  status_.num_jacobian_evaluations++;
  // must calculate derivatives with
  CalculateDTotal();

  // zero Jacobian
  jacobian_.Zero();
  // add in derivatives for equilibrium chemistry
  AddAccumulationDerivative(&jacobian_, &dtotal_, &dtotal_sorbed_);

  // add in derivatives for kinetic chemistry
  AddKineticChemistryToJacobian();
}


void Beaker::ScaleRHSAndJacobian() {
  for (int i = 0; i < jacobian_.size(); i++) {
    double max = jacobian_.GetRowAbsMax(i);
    if (max > 1.0) {
      double scale = 1.0 / max;
      rhs_.at(i) *= scale;
      jacobian_.ScaleRow(i, scale);
    }
  }
}


void Beaker::UpdateMolalitiesWithTruncation(const double max_ln_change) {
  double max_linear_change = std::pow(10.0, max_ln_change); // log10 vs ln... close enough
  double max_change;
  if (use_log_formulation_) {
    max_change = max_ln_change;
  } else {
    max_change = max_linear_change;
  }
  
  double min_ratio = 1.0e20; // large number

  for (int i = 0; i < ncomp_; i++) {
    // truncate the rhs to max_change
    if (rhs_.at(i) > max_change) {
      rhs_.at(i) = max_change;
    } else if (rhs_.at(i) < -max_change) {
      rhs_.at(i) = -max_change;
    }

    // store the previous solution
    prev_molal_.at(i) = primary_species().at(i).molality();

    if (!use_log_formulation_) {
      // ensure non-negative concentration
      if (prev_molal_.at(i) <= rhs_.at(i)) {
        double ratio = std::fabs(prev_molal_.at(i) / rhs_.at(i));
        min_ratio = std::min(ratio, min_ratio);
      }
    }
  }

  // update primary species molalities (log formulation)
  for (int i = 0; i < ncomp_; ++i) {
    double molality;
    if (use_log_formulation_) {
      molality = prev_molal_.at(i) * std::exp(-rhs_.at(i));
    } else {
      if (min_ratio < 1.0) {
        // scale by 0.99 to make the update slightly smaller than the min_ratio
        rhs_.at(i) *= min_ratio*0.99;
      }
      molality = prev_molal_.at(i) - rhs_.at(i);
    }
    primary_species_.at(i).update(molality);
  }
}


void Beaker::CalculateMaxRelChangeInMolality(double* max_rel_change, int* max_rel_index) {
  *max_rel_change = 0.0;
  *max_rel_index = -1;
  for (int i = 0; i < ncomp_; i++) {
    double delta = std::fabs(primary_species().at(i).molality() - prev_molal_.at(i)) / prev_molal_.at(i);
    if (delta > *max_rel_change) {
      *max_rel_change = delta;
      *max_rel_index = i;
    }
  }
} 


void Beaker::CheckChargeBalance(const std::vector<double>& aqueous_totals) const {
  double charge_balance = 0.0;
  for (unsigned int i = 0; i < aqueous_totals.size(); i++) {
    charge_balance += aqueous_totals.at(i) * primary_species().at(i).charge();
  }
  if (std::fabs(charge_balance) > tolerance_) {
    std::stringstream message;
    message << "WARNING: Beaker::CheckChargeBalance() : " 
            << " charge balance = " << std::scientific
            << charge_balance << std::fixed << std::endl;
    vo_->WriteWarning(Teuchos::VERB_EXTREME, message);
  }
}


void Beaker::ValidateSolution() {
  // TODO(bandre): what checks can we to to verify that the current solution is good?
  // TODO(bandre): check for nan or inf's as a sign that the step was too big?

  // TODO(bandre): negative total's (H+) are OK...

  // charge balance is error or warning...?
  CheckChargeBalance(total_);

  // negative mineral volume fractions are bad...
  for (unsigned int m = 0; m < minerals_.size(); m++) {
    if (minerals_.at(m).volume_fraction() < 0.0) {
      std::ostringstream error_stream;
      error_stream << "Beaker::ValidateSolution(): \n";
      error_stream << "   mineral " << minerals_.at(m).name()
                   << " volume_fraction is negative: " 
                   << minerals_.at(m).volume_fraction() << std::endl;
      Exceptions::amanzi_throw(ChemistryInvalidSolution(error_stream.str()));
    }
  }
}


/*******************************************************************************
 **
 **  Output related functions
 **
 *******************************************************************************/
void Beaker::DisplayParameters() const {
  std::stringstream message;
  // units....
  message << "---- Parameters" << std::endl;
  // message << "    thermo_database_file: " << thermo_database_file << std::endl;
  message << "    tolerance: " << tolerance_ << std::endl;
  message << "    max_iterations :" << max_iterations_ << std::endl;

  message << "    activity model: " << activity_model_->name() << std::endl;

  message << "    porosity: " << porosity_ << " [-]" << std::endl;
  message << "    water saturation: " << saturation_ << " [-]" << std::endl;
  message << "    water density: " << water_density_kg_m3() << " [kg m^-3]" << std::endl;
  message << "    volume: " << volume_ << " [m^3]" << std::endl;
  message << std::endl;
  vo_->Write(Teuchos::VERB_HIGH, message.str());
}


void Beaker::DisplayPrimary() const {
  std::stringstream message;
  message << "---- Primary Species" << std::endl;
  message << std::setw(15) << "Species"
          << std::setw(10) << "Charge"
          << std::setw(10) << "GMW"
          << std::setw(10) << "D-H a0"
          << std::endl;
  vo_->Write(Teuchos::VERB_HIGH, message.str());
  for (std::vector<Species>::const_iterator primary = primary_species().begin();
       primary != primary_species().end(); primary++) {
    primary->Display(vo_);
  }
  vo_->Write(Teuchos::VERB_HIGH, "\n");
}


void Beaker::DisplayAqueousEquilibriumComplexes() const {
  std::stringstream message;
  message << "---- Aqueous Equilibrium Complexes" << std::endl;
  message << std::setw(12) << "Reaction"
          << std::setw(38) << "log Keq"
          << std::setw(8) << "Charge"
          << std::setw(10) << "GMW"
          << std::setw(8) << "D-H a0"
          << std::endl;
  vo_->Write(Teuchos::VERB_HIGH, message.str());
  for (std::vector<AqueousEquilibriumComplex>::const_iterator aec = aqComplexRxns_.begin();
       aec != aqComplexRxns_.end(); aec++) {
    aec->Display(vo_);
  }
  vo_->Write(Teuchos::VERB_HIGH, "\n");
}


void Beaker::DisplayGeneralKinetics() const {
  if (generalKineticRxns_.size() > 0) {
    std::stringstream message;
    message << "---- General Kinetics" << std::endl;
    message << std::setw(12) << "Reaction" << std::endl;
    vo_->Write(Teuchos::VERB_HIGH, message.str());
    for (std::vector<GeneralRxn>::const_iterator rxn = generalKineticRxns_.begin();
         rxn != generalKineticRxns_.end(); rxn++) {
      rxn->Display(vo_);
    }
    vo_->Write(Teuchos::VERB_HIGH, "\n");
  }
}


void Beaker::DisplayRadioactiveDecayRxns() const {
  if (radioactive_decay_rxns_.size() > 0) {
    std::stringstream message;
    message << "---- Radioactive Decay" << std::endl;
    message << std::setw(12) << "Reaction" << std::endl;
    vo_->Write(Teuchos::VERB_HIGH, message.str());
    std::vector<RadioactiveDecay>::const_iterator rxn;
    for (rxn = radioactive_decay_rxns_.begin();
         rxn != radioactive_decay_rxns_.end(); ++rxn) {
      rxn->Display(vo_);
    }
    vo_->Write(Teuchos::VERB_HIGH, "\n");
  }
}


void Beaker::DisplayMinerals() const {
  if (minerals_.size() > 0) {
    std::stringstream message;
    message << "---- Minerals" << std::endl;
    message << std::setw(12) << "Reaction"
            << std::setw(38) << "log_Keq"
            << std::setw(13) << "molar volume"
            << std::setw(13) << "GMW"
            << std::setw(13) << "SSA"
            << std::setw(13) << "Vfrac"
            << std::endl;
    message << std::setw(12) << " "
            << std::setw(38) << " "
            << std::setw(13) << "[m^3/mol]"
            << std::setw(13) << "[g/mol]"
            << std::setw(13) << "[m^2/m^3 blk]"
            << std::setw(13) << "[-]"
            << std::endl;
    vo_->Write(Teuchos::VERB_HIGH, message.str());
    for (std::vector<Mineral>::const_iterator m = minerals_.begin();
         m != minerals_.end(); m++) {
      m->Display(vo_);
    }
    vo_->Write(Teuchos::VERB_HIGH, "\n");
  }
} 


void Beaker::DisplayMineralKinetics() const {
  if (mineral_rates_.size() > 0) {
    std::stringstream message;
    vo_->Write(Teuchos::VERB_HIGH, "---- Mineral Kinetics\n");
    for (std::vector<KineticRate*>::const_iterator m = mineral_rates_.begin();
         m != mineral_rates_.end(); m++) {
      (*m)->Display(vo_);
    }
    vo_->Write(Teuchos::VERB_HIGH, "\n");
  }
} 


void Beaker::DisplayIonExchangeSites() const {
  if (ion_exchange_rxns_.size() > 0) {
    std::stringstream message;
    message << "---- Ion Exchange Sites" << std::endl;
    message << std::setw(15) << "Species"
            << std::setw(20) << "Location"
            << std::setw(10) << "Charge"
            << std::setw(10) << "CEC"
            << std::endl;
    vo_->Write(Teuchos::VERB_HIGH, message.str());
    std::vector<IonExchangeRxn>::const_iterator rxn;
    for (rxn = ion_exchange_rxns_.begin();
         rxn != ion_exchange_rxns_.end(); rxn++) {
      rxn->site().Display(vo_);
    }
    vo_->Write(Teuchos::VERB_HIGH, "\n");
  }
}


void Beaker::DisplayIonExchangeComplexes() const {
  if (ion_exchange_rxns_.size() > 0) {
    std::stringstream message;
    message << "---- Ion Exchange Complexes" << std::endl;
    message << std::setw(12) << "Reaction"
            << std::setw(38) << "K"
            << std::endl;
    vo_->Write(Teuchos::VERB_HIGH, message.str());
    std::vector<IonExchangeRxn>::const_iterator ier;
    for (ier = ion_exchange_rxns_.begin();
         ier != ion_exchange_rxns_.end(); ier++) {
      ier->Display(vo_);
    }
    vo_->Write(Teuchos::VERB_HIGH, "\n");
  }
}


void Beaker::DisplaySurfaceSites() const {
  if (surfaceComplexationRxns_.size() > 0) {
    std::stringstream message;
    message << "---- Surface Sites" << std::endl;
    message << std::setw(15) << "Species"
            << std::setw(15) << "Site Density"
            << std::endl;
    message << std::setw(15) << " "
            << std::setw(15) << "[mol/m^3]"
            << std::endl;
    vo_->Write(Teuchos::VERB_HIGH, message.str());
    std::vector<SurfaceComplexationRxn>::const_iterator s;
    for (s = surfaceComplexationRxns_.begin();
         s != surfaceComplexationRxns_.end(); s++) {
      s->DisplaySite(vo_);
    }
    vo_->Write(Teuchos::VERB_HIGH, "\n");
  }
}


void Beaker::DisplaySurfaceComplexes() const {
  if (surfaceComplexationRxns_.size() > 0) {
    std::stringstream message;
    message << "---- Surface Complexes" << std::endl;
    message << std::setw(12) << "Reaction"
            << std::setw(38) << "log Keq"
            << std::setw(10) << "charge"
            << std::endl;
    vo_->Write(Teuchos::VERB_HIGH, message.str());
    std::vector<SurfaceComplexationRxn>::const_iterator s;
    for (s = surfaceComplexationRxns_.begin();
         s != surfaceComplexationRxns_.end(); s++) {
      s->DisplayComplexes(vo_);
    }
    vo_->Write(Teuchos::VERB_HIGH, "\n");
  }
}


void Beaker::DisplaySorptionIsotherms() const {
  if (sorption_isotherm_rxns_.size() > 0) {
    std::stringstream message;
    message << "---- Equilibrium Sorption Isotherms" << std::endl;
    message << std::setw(12) << "Species"
            << std::setw(15) << "isotherm"
            << std::setw(15) << "parameters"
            << std::endl;
    vo_->Write(Teuchos::VERB_HIGH, message.str());
    std::vector<SorptionIsothermRxn>::const_iterator s;
    for (s = sorption_isotherm_rxns_.begin();
         s != sorption_isotherm_rxns_.end(); s++) {
      s->Display(vo_);
    }
    vo_->Write(Teuchos::VERB_HIGH, "\n");
  }
}


void Beaker::print_results() const {
  // output for testing purposes
  std::stringstream message;
  message << std::endl;
  message << "----- Solution ----------------------" << std::endl;
  message << "Primary Species ---------------------\n";
  for (int i = 0; i < ncomp_; i++) {
    message << "   " << primary_species().at(i).name() << std::endl;
    message << "        Total: " << total_.at(i) << std::endl;
    message << "     Free-Ion: " << primary_species().at(i).molality() << std::endl;
    message << "Activity Coef: " << primary_species().at(i).act_coef() << std::endl;
    message << "     Activity: " << primary_species().at(i).activity() << std::endl;
  }
  message << std::endl;
  message << "Secondary Species -------------------\n";
  for (unsigned int i = 0; i < aqComplexRxns_.size(); i++) {
    message << "   " << aqComplexRxns_.at(i).name() << std::endl;
    message << "     Free-Ion: " << aqComplexRxns_.at(i).molality() << std::endl;
    message << "Activity Coef: " << aqComplexRxns_.at(i).act_coef() << std::endl;
    message << "     Activity: " << aqComplexRxns_.at(i).activity() << std::endl;
  }
  message << "-------------------------------------\n";
  message << std::endl;
}


void Beaker::print_results(double time) const {
  std::stringstream message;
  if (time < 1.e-40) {
    message << "Time\t";
    for (int i = 0; i < ncomp_; i++) {
      message << primary_species().at(i).name() << " (total)\t";
      message << primary_species().at(i).name() << " (free-ion)\t";
    }
    message << std::endl;
  }
  // output for testing purposes
  message << time << "\t";
  for (int i = 0; i < ncomp_; i++) {
    message << total_.at(i) << "\t";
    message << primary_species().at(i).molality() << "\t";
  }
  message << std::endl;
} 

}  // namespace AmanziChemistry
}  // namespace Amanzi
