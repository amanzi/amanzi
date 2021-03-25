/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Driver class for evaluating geochemical related processes at a
  single computational node

  TODO(bandre): update mineral volume fractions to components.minerals after kinetics...
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


/*******************************************************************************
 **
 ** public interface for drivers (batch_chem, unstructured and
 ** structured process kernels)
 **
 ******************************************************************************/

//
// public setup related
//

void Beaker::Setup(const Beaker::BeakerComponents& components,
                   const Beaker::BeakerParameters& parameters) {
  SetParameters(parameters);

  SetupActivityModel(parameters.activity_model_name,
                     parameters.pitzer_database, parameters.jfunction_pitzer);
  ResizeInternalMemory(static_cast<int>(primary_species().size()));
  VerifyComponentSizes(components);
}


Beaker::BeakerParameters Beaker::GetDefaultParameters(void) const {
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


Beaker::BeakerParameters Beaker::GetCurrentParameters(void) const {
  // Extract data from the various chemistry objects and store them
  // into a beaker parameters object that can be handed off to the
  // driver.
  Beaker::BeakerParameters parameters;

  parameters.thermo_database_file.clear();

  parameters.tolerance = tolerance();
  parameters.max_iterations = max_iterations();

  parameters.activity_model_name.clear();
  if (activity_model_ != NULL) {
    parameters.activity_model_name = activity_model_->name();
  }

  parameters.porosity = porosity();
  parameters.saturation = saturation();
  parameters.water_density = water_density_kg_m3();  // kg / m^3
  parameters.volume = volume();  // m^3

  return parameters;
}


void Beaker::SetParameters(const Beaker::BeakerParameters& parameters) {
  // Take a parameters object that was created by the driver, and map
  // the data into the appropriate chemistry object, potentially over
  // riding some of our internal database data.
  tolerance(parameters.tolerance);
  max_iterations(parameters.max_iterations);
  porosity(parameters.porosity);
  water_density_kg_m3(parameters.water_density);  // den = [kg/m^3]
  saturation(parameters.saturation);
  volume(parameters.volume);  // vol = [m^3]
  update_accumulation_coefficients();
  update_por_sat_den_vol();
}


//
// public "computation engine" routines
//

// if no water density provided, default is 1000.0 kg/m^3
int Beaker::Speciate(Beaker::BeakerComponents* components,
                     const Beaker::BeakerParameters& parameters) {
  std::stringstream message;
  double speciation_tolerance = 1.e-12;
  double residual_tolerance = 1.e-12;
  ResetStatus();
  UpdateParameters(parameters, 1.0);  // NOTE: need dt=1 to avoid divide by zero
  CheckChargeBalance(components->total);

  CopyComponentsToBeaker(*components);

  // store current molalities
  for (int i = 0; i < ncomp(); i++) {
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
    for (int i = 0; i < ncomp(); i++) {
      residual_.at(i) = total_.at(i) - components->total.at(i);
    }
    // add derivatives of total with respect to free to Jacobian
    // units of Jacobian: kg water/sec
    jacobian_.Zero();
    CalculateDTotal();
    jacobian_.AddValues(&dtotal_);

    for (int i = 0; i < ncomp(); i++) {
      rhs_.at(i) = residual_.at(i);
    }

    // scale the Jacobian
    ScaleRHSAndJacobian();

    // for derivatives with respect to ln concentration, scale columns
    // by primary species concentrations
    for (int i = 0; i < ncomp(); i++) {
      jacobian_.ScaleColumn(i, primary_species().at(i).molality());
    }

    // call solver
    lu_solver_.Solve(&jacobian_, &rhs_);

    // calculate update truncating at a maximum of 5 in log space
    UpdateMolalitiesWithTruncation(5.0);
    // calculate maximum relative change in concentration over all species
    CalculateMaxRelChangeInMolality(&max_rel_change, &max_rel_index);

    /*
      for (int i = 0; i < ncomp(); i++) {
        message << primary_species().at(i).name() << " "
                << primary_species().at(i).molality() << " " << total_.at(i) << "\n";
      }
      vo_->Write(Teuchos::VERB_EXTREME, "");
    */

    num_iterations++;

    max_residual = 0.;
    for (int i = 0; i < ncomp(); i++)
      if (std::fabs(residual_.at(i)) > max_residual) {
        max_residual = std::fabs(residual_.at(i));
      }

    // if max_rel_change small enough, turn on activity coefficients
    if (max_rel_change < speciation_tolerance ||
        max_residual < residual_tolerance) {
      calculate_activity_coefs = true;
    }

    // exist if maximum relative change is below tolerance
  } while (max_rel_change > speciation_tolerance &&
           num_iterations < max_iterations());
  
  // for now, initialize total sorbed concentrations based on the current free
  // ion concentrations
  UpdateEquilibriumChemistry();
  CopyBeakerToComponents(components);
  status_.num_newton_iterations = num_iterations;
  if (max_rel_change < tolerance()) {
    status_.converged = true;
  }

  /*
    message.str("");
    message << "Beaker::speciate: max_rel_change: " << max_rel_change 
            << "  tolerance: " << speciation_tolerance << std::endl;
    message << "Beaker::speciate: max_residual: " << max_residual
            << "  tolerance: " << residual_tolerance << std::endl;
    message << "Beaker::speciate: status.num_rhs_evaluations: " << status_.num_rhs_evaluations << std::endl;
    message << "Beaker::speciate: status.num_jacobian_evaluations: " << status_.num_jacobian_evaluations << std::endl;
    message << "Beaker::speciate: status.num_newton_iterations: " << status_.num_newton_iterations << std::endl;
    message << "Beaker::speciate: status.converged: " << status_.converged << std::endl;
    vo_->Write(Teuchos::VERB_HIGH, message);
  */

  return num_iterations;
}


int Beaker::ReactionStep(Beaker::BeakerComponents* components,
                         const Beaker::BeakerParameters& parameters,
                         double dt) {
  /*
  ** Note: the parameter components is modified by this function.
  ** initially it contains the initial component concentrations.
  ** on return it contains the modified values of the components.
  */
  std::stringstream message;
  // update class paramters
  // water_density [kg/m^3]
  // volume [m^3]
  ResetStatus();
  UpdateParameters(parameters, dt);
  CheckChargeBalance(components->total);
  CopyComponentsToBeaker(*components);

  // store current molalities
  for (int i = 0; i < ncomp(); i++) {
    prev_molal_.at(i) = primary_species().at(i).molality();
  }

  // initialize to a large number (not necessary, but safe)
  double max_rel_change = 1.e20;
  int max_rel_index = -1;
  // iteration counter
  unsigned int num_iterations = 0;

  //set_use_log_formulation(false);

  // lagging activity coefficients by a time step in this case
  //UpdateActivityCoefficients();

  // calculate portion of residual at time level t
  CalculateFixedAccumulation(components->total, components->total_sorbed,
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

    for (int i = 0; i < ncomp(); i++) {
      rhs_.at(i) = residual_.at(i);
    }

    // scale the Jacobian
    ScaleRHSAndJacobian();

    if (use_log_formulation()) {
      // for derivatives with respect to ln concentration, scale columns
      // by primary species concentrations  
      for (int i = 0; i < ncomp(); i++) {
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

    /*
      message.str("");
      message << "--Iteration: "
              << num_iterations << "\n  max_rel_change(" << max_rel_index 
              << ") : " << max_rel_change << "\n";
      message << std::scientific << std::setprecision(16);
      for (int i = 0; i < ncomp(); i++) {
        message << std::setw(10) << primary_species().at(i).name()
                << std::setw(25) << primary_species().at(i).molality()
                << std::setw(25) << total_.at(i);
        if (total_sorbed_.size() > 0) {
          message << std::setw(25) << total_sorbed_.at(i);
        }
        message << "\n";
      }
      vo_->Write(Teuchos::VERB_HIGH, message);
    */

    num_iterations++;

    // geh
    //    if (num_iterations >= 100) {
    //      for (int i = 0; i < ncomp(); i++)
    //        std::cout << primary_species_.at(i).name() << " " <<
    //                  primary_species_.at(i).molality() << " " << total_.at(i) << "\n";
    //      std::cout << max_rel_change << " " << tolerance() << std::endl;
    //    }

    // exit if maximum relative change is below tolerance
  } while (max_rel_change > tolerance() && num_iterations < max_iterations());

  if (num_iterations >= max_iterations()) {
    // TODO(bandre): should this be an error to the driver...?
    // code eventually produces nans when this isn't an error.
    std::ostringstream error_stream;
    error_stream << "Warning: The maximum number Netwon iterations reached in Beaker::ReactionStep()." << std::endl;
    error_stream << "Warning: Results may not have the desired accuracy." << std::endl;
    error_stream << "Warning: max relative change = " << max_rel_change << std::endl;
    error_stream << "Warning: max relative index = " << max_rel_index << std::endl;
    error_stream << "Warning: tolerance = " << tolerance() << std::endl;
    error_stream << "Warning: max iterations = " << max_iterations() << std::endl;
    // update before leaving so that we can see the erroneous values!
    CopyBeakerToComponents(components);
    Exceptions::amanzi_throw(ChemistryMaxIterationsReached(error_stream.str()));
  }

  status_.num_newton_iterations = num_iterations;
  if (max_rel_change < tolerance()) {
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

  CopyBeakerToComponents(components);
  ValidateSolution();

  return num_iterations;
}


void Beaker::CopyBeakerToComponents(Beaker::BeakerComponents* components) {
  // Copy the beaker state into variables are be returned to the
  // driver for long term storage.

  // NOTE: The components struct may have been only partially
  // initialized be the driver (it doesn't know about the size of
  // internal memory for surface complexation or ion exchange....
  // For now we assert the things the driver should know....

  //
  // totals
  //
  assert(components->total.size() == total_.size());
  for (int i = 0; i < ncomp(); ++i) {
    components->total.at(i) = total_.at(i);
  }

  if (total_sorbed_.size() > 0) {
    assert(components->total_sorbed.size() == total_sorbed_.size());
    for (int i = 0; i < ncomp(); ++i) {
      components->total_sorbed.at(i) = total_sorbed_.at(i);
    }
  }

  //
  // free ion
  //
  assert(components->free_ion.size() == ncomp());
  for (int i = 0; i < ncomp(); ++i) {
    components->free_ion.at(i) = primary_species().at(i).molality();
  }

  //
  // activity coeff
  //
  if (components->primary_activity_coeff.size() != ncomp()) {
    components->primary_activity_coeff.resize(ncomp());
  }
  for (int i = 0; i < ncomp(); ++i) {
    components->primary_activity_coeff.at(i) = primary_species().at(i).act_coef();
  }
  if (components->secondary_activity_coeff.size() != aqComplexRxns_.size()) {
    components->secondary_activity_coeff.resize(aqComplexRxns_.size());
  }
  for (int i = 0; i < aqComplexRxns_.size(); ++i) {
    components->secondary_activity_coeff.at(i) = aqComplexRxns_.at(i).act_coef();
  }

  //
  // minerals
  //
  assert(components->mineral_volume_fraction.size() == minerals_.size());
  for (unsigned int m = 0; m < minerals_.size(); ++m) {
    components->mineral_volume_fraction.at(m) = minerals_.at(m).volume_fraction();
  }

  if (components->mineral_specific_surface_area.size() != minerals_.size()) {
    components->mineral_specific_surface_area.resize(minerals_.size());
  }
  for (int m = 0; m < minerals_.size(); ++m) {
    double ssa = minerals_.at(m).specific_surface_area();
    components->mineral_specific_surface_area.at(m) = ssa;
  }

  //
  // ion exchange
  //
  if (components->ion_exchange_sites.size() != ion_exchange_rxns_.size()) {
    components->ion_exchange_sites.resize(ion_exchange_rxns_.size());
  }
  for (int i = 0; i < ion_exchange_rxns_.size(); ++i) {
    components->ion_exchange_sites.at(i) = 
        ion_exchange_rxns_.at(i).site().get_cation_exchange_capacity();
  }

  if (components->ion_exchange_ref_cation_conc.size() != ion_exchange_rxns_.size()) {
    components->ion_exchange_ref_cation_conc.resize( 
        ion_exchange_rxns_.size());
  }
  for (int i = 0; i < ion_exchange_rxns_.size(); ++i) {
    components->ion_exchange_ref_cation_conc.at(i) = 
        ion_exchange_rxns_[i].ref_cation_sorbed_conc();
  }

  //
  // surface complexation
  //
  if (components->surface_complex_free_site_conc.size() != 
      surfaceComplexationRxns_.size()) {
    components->surface_complex_free_site_conc.resize( 
        surfaceComplexationRxns_.size());
  }
  for (unsigned int i = 0; i < surfaceComplexationRxns_.size(); ++i) {
    components->surface_complex_free_site_conc.at(i) =
        surfaceComplexationRxns_.at(i).free_site_concentration();
  }

  if (components->surface_site_density.size() != 
      surfaceComplexationRxns_.size()) {
    components->surface_site_density.resize(surfaceComplexationRxns_.size());
  }
  for (unsigned int i = 0; i < surfaceComplexationRxns_.size(); ++i) {
    components->surface_site_density.at(i) =
        surfaceComplexationRxns_.at(i).GetSiteDensity();
  }

  //
  // sorption isotherms
  //
  if (sorption_isotherm_rxns_.size() > 0) {
    if (components->isotherm_kd.size() != ncomp()) {
      components->isotherm_kd.resize(ncomp(), 0.0);
    }
    if (components->isotherm_langmuir_b.size() != ncomp()) {
      components->isotherm_langmuir_b.resize(ncomp(), 0.0);
    }
    if (components->isotherm_freundlich_n.size() != ncomp()) {
      components->isotherm_freundlich_n.resize(ncomp(), 1.0);
    }
    for (int r = 0; r < sorption_isotherm_rxns_.size(); ++r) {
      const std::vector<double>& params = sorption_isotherm_rxns_.at(r).GetIsothermParameters();
      int id = sorption_isotherm_rxns_.at(r).species_id();
      components->isotherm_kd.at(id) = params.at(0);
      if (sorption_isotherm_rxns_.at(r).IsothermType() == SorptionIsotherm::FREUNDLICH) {
        components->isotherm_freundlich_n.at(id) = params.at(1);
      } else if (sorption_isotherm_rxns_.at(r).IsothermType() == SorptionIsotherm::LANGMUIR) {
        components->isotherm_langmuir_b.at(id) = params.at(1);
      }
    }
  }
}


void Beaker::CopyComponents(const Beaker::BeakerComponents& from,
                            Beaker::BeakerComponents* to) {
  // this function doesn't really do anything...
  *to = from;
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

void Beaker::Display(void) const {
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


void Beaker::DisplayComponents(const Beaker::BeakerComponents& components) const {
  std::stringstream message;
  message << "--- Input Components -------------------------------------------------"
          << std::endl;
  message << "---- Aqueous Components" << std::endl;
  message << std::setw(15) << "Name"
          << std::setw(15) << "Molality"
          << std::setw(15) << "Molarity"
          // << std::setw(15) << "Free Ion" // TODO(bandre): uncomment and update test results
          << std::endl;
  for (int i = 0; i < ncomp(); i++) {
    message << std::setw(15) << primary_species().at(i).name()
            << std::scientific << std::setprecision(5)
            << std::setw(15) << components.total.at(i) / water_density_kg_L()
            << std::setw(15) << components.total.at(i)
            // << std::setw(15) << components.free_ion.at(i) // TODO(bandre): uncomment and update test results
            << std::endl;
  }

  if (minerals_.size() > 0) {
    message << "---- Mineral Components" << std::endl;
    message << std::setw(15) << "Name"
            << std::setw(15) << "Vol. frac" << std::endl;
    for (unsigned int m = 0; m < minerals_.size(); m++) {
      message << std::setw(15) << minerals_.at(m).name()
              << std::setw(15) << std::fixed << std::setprecision(5)
              << components.mineral_volume_fraction.at(m) << std::endl;
    }
  }

  if (total_sorbed_.size() > 0) {
    message << "---- Sorbed Components" << std::endl;
    message << std::setw(15) << "Name"
            << std::setw(15) << "Moles / m^3" << std::endl;
    for (int i = 0; i < ncomp(); i++) {
      message << std::setw(15) << primary_species().at(i).name()
              << std::scientific << std::setprecision(5)
              << std::setw(15) << components.total_sorbed.at(i)
              << std::endl;
    }
  }
  message << "------------------------------------------------- Input Components ---"
          << std::endl;
  vo_->Write(Teuchos::VERB_HIGH, message.str());
}


void Beaker::DisplayResults(void) const {
  std::stringstream message;
  message << std::endl;
  message << "-- Solution ----------------------------------------------------------"
          << std::endl;
  message << "---- Components " << std::endl;
  message << std::setw(15) << "Name"
          << std::setw(15) << "Molality"
          << std::setw(15) << "Molarity"
          << std::endl;
  for (int i = 0; i < ncomp(); i++) {
    message << std::setw(15) << primary_species().at(i).name()
            << std::scientific << std::setprecision(5)
            << std::setw(15) << total_.at(i) / water_density_kg_L()
            << std::setw(15) << total_.at(i)
            << std::endl;
  }

  message << "---- Change Balance " << std::endl;
  double charge_balance_molal = 0.0;
  for (int i = 0; i < ncomp(); i++) {
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
  for (int i = 0; i < ncomp(); i++) {
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
  for (int i = 0; i < ncomp(); i++) {
    message << std::setw(15) << primary_species().at(i).name();
  }
  if (display_free) {
    for (int i = 0; i < ncomp(); i++) {
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
                                 const BeakerComponents& components,
                                 const bool display_free) const {
  std::stringstream message;
  message << std::scientific << std::setprecision(6) << std::setw(15);
  message << time;
  for (int i = 0; i < ncomp(); i++) {
    message << std::setw(15) << components.total.at(i);
  }
  if (display_free) {
    for (int i = 0; i < ncomp(); i++) {
      message << std::setw(15) << components.free_ion.at(i);
    }
  }
  if (total_sorbed_.size() > 0) {
    for (int i = 0; i < total_sorbed_.size(); i++) {
      message << std::setw(15) << components.total_sorbed.at(i);
    }
  }
  if (minerals().size() > 0) {
    for (int m = 0; m < minerals().size(); ++m) {
      message << std::setw(15) << components.mineral_volume_fraction.at(m);
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
  ncomp(size);
  total_.resize(ncomp());
  dtotal_.Resize(ncomp());

  if (surfaceComplexationRxns_.size() > 0 ||
      sorption_isotherm_rxns_.size() > 0 ||
      ion_exchange_rxns_.size() > 0) {
    total_sorbed_.resize(ncomp(), 0.0);
    dtotal_sorbed_.Resize(ncomp());
    dtotal_sorbed_.Zero();
  } else {
    total_sorbed_.resize(0);
    // dtotal_sorbed_.Resize(0);
  }

  fixed_accumulation_.resize(ncomp());
  residual_.resize(ncomp());
  prev_molal_.resize(ncomp());

  jacobian_.Resize(ncomp());
  rhs_.resize(ncomp());
  lu_solver_.Initialize(ncomp());
}


void Beaker::VerifyComponentSizes(const Beaker::BeakerComponents& components) const {
  // some helpful error checking goes here...
  bool error = false;
  std::ostringstream error_stream;
  error_stream << "Beaker::VerifyComponentSizes():\ndatabase input and component initial conditions do not match:\n";

  // if the size of the various initial conditions, components, and
  // database input don't match. Print a helpful message and exit
  // gracefully.

  if (static_cast<unsigned int>(this->ncomp()) != components.total.size()) {
    error = true;
    error_stream << "ncomp(" << this->ncomp()
                 << ") and components.total.size(" << components.total.size()
                 << ") do not match.\n";
  }

  if (this->primary_species().size() != components.total.size()) {
    error = true;
    error_stream << "primary_species.size(" << this->primary_species().size()
                 << ") and components.total.size(" << components.total.size()
                 << ") do not match.\n";
  }

  if (components.free_ion.size() != components.total.size()) {
    error = true;
    error_stream << "components.total.size(" << components.total.size()
                 << ") and components.free_ion.size(" << components.free_ion.size()
                 << ") do not match.\n";
  }

  if (this->ion_exchange_rxns().size() != components.ion_exchange_sites.size()) {
    error = true;
    error_stream << "ion_exchange_rxns.size(" << this->ion_exchange_rxns().size()
                 << ") and components.ion_exchange_sites.size(" << components.ion_exchange_sites.size()
                 << ") do not match.\n";
  }

  if (this->minerals().size() != components.mineral_volume_fraction.size()) {
    error = true;
    error_stream << "minerals.size(" << this->minerals().size()
                 << ") and components.mineral_volume_fraction.size(" << components.mineral_volume_fraction.size()
                 << ") do not match.\n";
  }

  if (this->total().size() != components.total.size()) {
    error = true;
    error_stream << "total.size(" << this->total().size()
                 << ") and components.total.size(" << components.total.size()
                 << ") do not match.\n";
  }

  // FIXED? this check is breaking things because total_sorbed is always
  // resized in resize(), even if there is no sorption!
  if (this->total_sorbed().size() != components.total_sorbed.size()) {
    error = true;
    error_stream << "total_sorbed.size(" << this->total_sorbed().size()
                 << ") and components.total_sorbed.size("
                 << components.total_sorbed.size() << ") do not match.\n";
  }

  if (error) {
    Exceptions::amanzi_throw(ChemistryMemorySizeError(error_stream.str()));
  }
}


void Beaker::CopyComponentsToBeaker(const Beaker::BeakerComponents& components) {
  // NOTE: Do not copy total and total_sorbed here!

  //
  // free ion
  //
  if (components.free_ion.size() > 0) {
    InitializeMolalities(components.free_ion);
  } else {
    InitializeMolalities(1.e-9);
  }

  //
  // activity coefficients
  //
  if (components.primary_activity_coeff.size() > 0) {
    assert(components.primary_activity_coeff.size() == primary_species().size());
    for (int i = 0; i < primary_species().size(); ++i) {
      double value = components.primary_activity_coeff.at(i);
      primary_species_.at(i).act_coef(value);
    }
  }
  if (components.secondary_activity_coeff.size() > 0) {
    assert(components.secondary_activity_coeff.size() == aqComplexRxns_.size());
    for (int i = 0; i < aqComplexRxns_.size(); ++i) {
      double value = components.secondary_activity_coeff.at(i);
      aqComplexRxns_.at(i).act_coef(value);
    }
  }

  //
  // minerals
  //
  unsigned int size = components.mineral_volume_fraction.size();
  if (minerals().size() == size) {
    for (unsigned int m = 0; m < size; m++) {
      minerals_.at(m).set_volume_fraction(components.mineral_volume_fraction.at(m));
    }
  } else {
    std::ostringstream error_stream;
    error_stream << "Beaker::CopyComponentsToBeaker(): \n";
    error_stream << "minerals.size(" << minerals().size()
                 << ") and components.mineral_volume_fraction.size(" << size
                 << ") do not match.\n";
    Exceptions::amanzi_throw(ChemistryUnrecoverableError(error_stream.str()));
  }

  if (components.mineral_specific_surface_area.size() > 0) {
    assert(components.mineral_specific_surface_area.size() == minerals_.size());
    for (unsigned int m = 0; m < minerals_.size(); ++m) {
      double value = components.mineral_specific_surface_area.at(m);
      minerals_.at(m).set_specific_surface_area(value);
    }
  }

  //
  // ion exchange
  //
  size = components.ion_exchange_sites.size();
  if (ion_exchange_rxns().size() == size) {
    for (unsigned int ies = 0; ies < size; ies++) {
      ion_exchange_rxns_[ies].set_cation_exchange_capacity(components.ion_exchange_sites.at(ies));
    }
  } else {
    std::ostringstream error_stream;
    error_stream << "Beaker::CopyComponentsToBeaker(): \n";
    error_stream << "ion_exchange_sites.size(" << ion_exchange_rxns().size()
                 << ") and components.ion_exchange_sites.size(" << size 
                 << ") do not match.\n";
    Exceptions::amanzi_throw(ChemistryUnrecoverableError(error_stream.str()));
  }

  if (ion_exchange_rxns().size() > 0) {
    // check for the ref cation conc....
    if (components.ion_exchange_ref_cation_conc.size() ==
        ion_exchange_rxns().size()) {
      // we have a value from a previous solve, restore it.
      for (int i = 0; i < ion_exchange_rxns().size(); ++i) {
        double value = components.ion_exchange_ref_cation_conc.at(i);
        ion_exchange_rxns_[i].set_ref_cation_sorbed_conc(value);
      }
    } else {
      // no previous value, provide a guess
      for (int i = 0; i < ion_exchange_rxns().size(); ++i) {
        ion_exchange_rxns_[i].set_ref_cation_sorbed_conc(1.0e-9);
      }
    }
  }

  //
  // surface complexation
  //
  if (components.surface_site_density.size() > 0) {
    assert(components.surface_site_density.size() == surfaceComplexationRxns_.size());
    for (unsigned int i = 0; i < surfaceComplexationRxns_.size(); ++i) {
      double value = components.surface_site_density.at(i);
      surfaceComplexationRxns_.at(i).UpdateSiteDensity(value);
    }
  }

  if (surfaceComplexationRxns_.size() > 0) {
    if (components.surface_complex_free_site_conc.size() ==
        surfaceComplexationRxns_.size()) {
      // we have a value from a previous solve, restore it.
      for (int r = 0; r < surfaceComplexationRxns_.size(); ++r) {
        double value = components.surface_complex_free_site_conc.at(r);
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

  //
  // sorption isotherms
  //
  if (sorption_isotherm_rxns_.size() > 0 &&
      components.isotherm_kd.size() > 0) {
    // the driver maybe attempting to over the database values, or the
    // components were resized by a call to CopyBeakerToComponents()
    assert(components.isotherm_kd.size() == ncomp());
    assert(components.isotherm_freundlich_n.size() == ncomp());
    assert(components.isotherm_langmuir_b.size() == ncomp());
    // NOTE(bandre): sorption_isotherm_params_ hard coded size=4,
    // current max parameters is 2
    for (int r = 0; r < sorption_isotherm_rxns_.size(); ++r) {
      int id = sorption_isotherm_rxns_.at(r).species_id();
      sorption_isotherm_params_.at(0) = components.isotherm_kd.at(id);
      SorptionIsotherm::SorptionIsothermType isotherm_type = sorption_isotherm_rxns_.at(r).IsothermType();
      if (isotherm_type == SorptionIsotherm::FREUNDLICH) {
        sorption_isotherm_params_.at(1) = components.isotherm_freundlich_n.at(id);
      } else if (isotherm_type == SorptionIsotherm::LANGMUIR) {
        sorption_isotherm_params_.at(1) = components.isotherm_langmuir_b.at(id);
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


bool Beaker::HaveKinetics(void) const {
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
  dt(delta_t);  // delta time = [sec]
  SetParameters(parameters);
} 


void Beaker::ResetStatus(void) {
  status_.num_rhs_evaluations = 0;
  status_.num_jacobian_evaluations = 0;
  status_.num_newton_iterations = 0;
  status_.converged = false;
}


void Beaker::update_accumulation_coefficients(void) {
  aqueous_accumulation_coef(porosity() * saturation() * volume() * 1000.0 / dt());
  sorbed_accumulation_coef(volume() / dt());
} 


void Beaker::update_por_sat_den_vol(void) {
  por_sat_den_vol(porosity() * saturation() * water_density_kg_m3() * volume());
}


void Beaker::UpdateActivityCoefficients(void) {
  activity_model_->CalculateIonicStrength(primary_species(),
                                          aqComplexRxns_);
  activity_model_->CalculateActivityCoefficients(&primary_species_,
                                                 &aqComplexRxns_,
                                                 &water_);
  for (std::vector<Species>::iterator i = primary_species_.begin();
       i != primary_species_.end(); i++) {
    i->update();
  }
}


void Beaker::UpdateKineticMinerals(void) {

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
    minerals_.at(i).UpdateVolumeFraction(kinetic_rate, dt());
    minerals_.at(i).UpdateSpecificSurfaceArea();
  }
}


void Beaker::InitializeMolalities(double initial_molality) {
  for (std::vector<Species>::iterator i = primary_species_.begin();
       i != primary_species_.end(); i++) {
    i->update(initial_molality);
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


void Beaker::UpdateEquilibriumChemistry(void) {
  //    calculateActivityCoefficients(-1);

  // update primary species activities
  for (std::vector<Species>::iterator primary = primary_species_.begin();
       primary != primary_species_.end(); primary++) {
    primary->update();
  }
  // calculated seconday aqueous complex concentrations
  for (std::vector<AqueousEquilibriumComplex>::iterator aqcplx =
           aqComplexRxns_.begin();
       aqcplx != aqComplexRxns_.end(); aqcplx++) {
    aqcplx->Update(primary_species(), water_);
  }

  // calculate mineral saturation states
  for (std::vector<Mineral>::iterator m = minerals_.begin();
       m != minerals_.end(); m++) {
    m->Update(primary_species(), water_);
  }
  // surface complexation
  for (std::vector<SurfaceComplexationRxn>::iterator srfcplx =
           surfaceComplexationRxns_.begin();
       srfcplx != surfaceComplexationRxns_.end(); srfcplx++) {
    srfcplx->Update(primary_species());
  }
  // sorption isotherms
  for (std::vector<SorptionIsothermRxn>::iterator i =
           sorption_isotherm_rxns_.begin();
       i != sorption_isotherm_rxns_.end(); i++) {
    i->Update(primary_species());
  }
  // add equilibrium ion exchange here?
  for (std::vector<IonExchangeRxn>::iterator ier = ion_exchange_rxns_.begin();
       ier != ion_exchange_rxns_.end(); ier++) {
    ier->Update(primary_species());
  }

  // calculate total component concentrations
  CalculateTotal();
}


void Beaker::CalculateTotal(void) {
  // add in primaries
  for (unsigned int i = 0; i < total_.size(); i++) {
    total_.at(i) = primary_species().at(i).molality();
  }

  // add in aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::iterator i = aqComplexRxns_.begin();
       i != aqComplexRxns_.end(); i++) {
    i->AddContributionToTotal(&total_);
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
  for (std::vector<SurfaceComplexationRxn>::iterator i =
           surfaceComplexationRxns_.begin();
       i != surfaceComplexationRxns_.end(); i++) {
    i->AddContributionToTotal(&total_sorbed_);
  }
  // add in isotherm contributions
  for (std::vector<SorptionIsothermRxn>::iterator i =
           sorption_isotherm_rxns_.begin();
       i != sorption_isotherm_rxns_.end(); i++) {
    i->AddContributionToTotal(&total_sorbed_);
  }
  // add ion exchange
  for (std::vector<IonExchangeRxn>::iterator ier = ion_exchange_rxns_.begin();
       ier != ion_exchange_rxns_.end(); ier++) {
    ier->AddContributionToTotal(&total_sorbed_);
  }
}


void Beaker::CalculateDTotal(void) {
  dtotal_.Zero();
  // derivative with respect to free-ion is 1.
  dtotal_.SetDiagonal(1.0);

  // add in derviative of complex contribution with respect to free-ion
  for (std::vector<AqueousEquilibriumComplex>::iterator i = aqComplexRxns_.begin();
       i != aqComplexRxns_.end(); i++) {
    i->AddContributionToDTotal(primary_species(), &dtotal_);
  }

  // scale by density of water
  dtotal_.Scale(water_density_kg_L());
  // dtotal_.Print("-- dtotal_ scaled");

  // calculate sorbed derivatives
  if (total_sorbed_.size()) {
    dtotal_sorbed_.Zero();
    for (std::vector<SurfaceComplexationRxn>::iterator i =
             surfaceComplexationRxns_.begin();
         i != surfaceComplexationRxns_.end(); i++) {
      i->AddContributionToDTotal(primary_species(), &dtotal_sorbed_);
    }
    for (std::vector<SorptionIsothermRxn>::iterator i =
             sorption_isotherm_rxns_.begin();
         i != sorption_isotherm_rxns_.end(); i++) {
      i->AddContributionToDTotal(primary_species(), &dtotal_sorbed_);
    }
    // add ion exchange
    for (std::vector<IonExchangeRxn>::iterator ier = ion_exchange_rxns_.begin();
         ier != ion_exchange_rxns_.end(); ier++) {
      ier->AddContributionToDTotal(primary_species(), &dtotal_sorbed_);
    }
  }
}


void Beaker::UpdateKineticChemistry(void) {
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
    i->UpdateRate(total_, total_sorbed_, porosity(), saturation(), volume());
  }

  // add mineral saturation and rate calculations here
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->Update(primary_species(), minerals_);
  }
  // add multirate kinetic surface complexation reaction quotient calculations
  // here
} 


void Beaker::AddKineticChemistryToResidual(void) {
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
    (*rate)->AddContributionToResidual(minerals_, volume(), &residual_);
  }

  // add multirate kinetic surface complexation contribution to residual here.
}


void Beaker::AddKineticChemistryToJacobian(void) {
  // loop over general kinetic reactions and add rates
  for (std::vector<GeneralRxn>::iterator i = generalKineticRxns_.begin();
       i != generalKineticRxns_.end(); i++) {
    i->addContributionToJacobian(&jacobian_, primary_species(), por_sat_den_vol());
  }

  // loop over radioactive decay reactions and add rates
  for (std::vector<RadioactiveDecay>::iterator i = radioactive_decay_rxns_.begin();
       i != radioactive_decay_rxns_.end(); ++i) {
    i->AddContributionToJacobian(dtotal_, dtotal_sorbed_,
                                 porosity(), saturation(), volume(),
                                 &jacobian_);
  }

  // add mineral mineral contribution to Jacobian here.  units = kg water/sec.
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->AddContributionToJacobian(primary_species(), minerals_, volume(), &jacobian_);
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
    residual->at(i) += aqueous_accumulation_coef() * total.at(i);
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
  J->AddValues(dtotal, aqueous_accumulation_coef());

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


void Beaker::CalculateResidual(void) {
  status_.num_rhs_evaluations++;
  // subtract fixed porition
  for (int i = 0; i < ncomp(); i++) {
    residual_.at(i) = -fixed_accumulation().at(i);
  }

  // accumulation adds in equilibrium chemistry
  AddAccumulation(total_, total_sorbed_, &residual_);

  // kinetic reaction contribution to residual
  AddKineticChemistryToResidual();
}


void Beaker::CalculateJacobian(void) {
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


void Beaker::ScaleRHSAndJacobian(void) {
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
  if (use_log_formulation()) {
    max_change = max_ln_change;
  } else {
    max_change = max_linear_change;
  }
  
  double min_ratio = 1.0e20; // large number

  for (int i = 0; i < ncomp(); i++) {
    // truncate the rhs to max_change
    if (rhs_.at(i) > max_change) {
      rhs_.at(i) = max_change;
    } else if (rhs_.at(i) < -max_change) {
      rhs_.at(i) = -max_change;
    }

    // store the previous solution
    prev_molal_.at(i) = primary_species().at(i).molality();

    if (!use_log_formulation()) {
      // ensure non-negative concentration
      if (prev_molal_.at(i) <= rhs_.at(i)) {
        double ratio = std::fabs(prev_molal_.at(i) / rhs_.at(i));
        min_ratio = std::min(ratio, min_ratio);
      }
    }  // if (use_log_formulation)
  }  // for (i)

  // update primary species molalities (log formulation)
  for (int i = 0; i < ncomp(); ++i) {
    double molality;
    if (use_log_formulation()) {
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
  for (int i = 0; i < ncomp(); i++) {
    double delta = std::fabs(primary_species().at(i).molality() - prev_molal().at(i)) / 
        prev_molal().at(i);
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
  if (std::fabs(charge_balance) > tolerance()) {
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
void Beaker::DisplayParameters(void) const {
  std::stringstream message;
  // units....
  message << "---- Parameters" << std::endl;
  // message << "    thermo_database_file: " << thermo_database_file << std::endl;
  message << "    tolerance: " << tolerance() << std::endl;
  message << "    max_iterations :" << max_iterations() << std::endl;

  message << "    activity model: " << activity_model_->name() << std::endl;

  message << "    porosity: " << porosity() << " [-]" << std::endl;
  message << "    water saturation: " << saturation() << " [-]" << std::endl;
  message << "    water density: " << water_density_kg_m3() << " [kg m^-3]" << std::endl;
  message << "    volume: " << volume() << " [m^3]" << std::endl;
  message << std::endl;
  vo_->Write(Teuchos::VERB_HIGH, message.str());
}


void Beaker::DisplayPrimary(void) const {
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


void Beaker::DisplayAqueousEquilibriumComplexes(void) const {
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


void Beaker::DisplayGeneralKinetics(void) const {
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


void Beaker::DisplayRadioactiveDecayRxns(void) const {
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


void Beaker::DisplayMinerals(void) const {
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


void Beaker::DisplayMineralKinetics(void) const {
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


void Beaker::DisplayIonExchangeSites(void) const {
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


void Beaker::DisplayIonExchangeComplexes(void) const {
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


void Beaker::DisplaySurfaceSites(void) const {
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


void Beaker::DisplaySurfaceComplexes(void) const {
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


void Beaker::DisplaySorptionIsotherms(void) const {
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


void Beaker::print_results(void) const {
  // output for testing purposes
  std::stringstream message;
  message << std::endl;
  message << "----- Solution ----------------------" << std::endl;
  message << "Primary Species ---------------------\n";
  for (int i = 0; i < ncomp(); i++) {
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
    for (int i = 0; i < ncomp(); i++) {
      message << primary_species().at(i).name() << " (total)\t";
      message << primary_species().at(i).name() << " (free-ion)\t";
    }
    message << std::endl;
  }
  // output for testing purposes
  message << time << "\t";
  for (int i = 0; i < ncomp(); i++) {
    message << total_.at(i) << "\t";
    message << primary_species().at(i).molality() << "\t";
  }
  message << std::endl;
} 

}  // namespace AmanziChemistry
}  // namespace Amanzi
