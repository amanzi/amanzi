/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** TODO(bandre): update mineral volume fractions to components.minerals after kinetics....
**
** TODO(bandre): finish implementing ion exchange jacobian
*/

#include "beaker.hh"

#include <cstdlib>
#include <cassert>

#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>

#include "activity_model.hh"
#include "activity_model_factory.hh"
#include "aqueous_equilibrium_complex.hh"
#include "general_rxn.hh"
#include "ion_exchange_rxn.hh"
#include "kinetic_rate.hh"
#include "mineral.hh"
#include "mineral_kinetics_factory.hh"
#include "species.hh"
#include "surface_complexation_rxn.hh"
#include "lu_solver.hh"
#include "matrix_block.hh"
#include "chemistry_verbosity.hh"
#include "chemistry_exception.hh"

#include "exceptions.hh"

namespace amanzi {
namespace chemistry {

extern ChemistryOutput* chem_out;

// solver defaults
const double Beaker::tolerance_default_ = 1.0e-12;
const unsigned int Beaker::max_iterations_default_ = 250;
// default physical parameters
const double Beaker::porosity_default_ = 1.0;  // [-]
const double Beaker::saturation_default_ = 1.0;  // [-]
const double Beaker::water_density_kg_m3_default_ = 1000.0;
const double Beaker::volume_default_ = 1.0;  // [m^3]

Beaker::Beaker()
    : debug_(false),
      verbosity_(kSilent),
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
      dt_(0.0),
      aqueous_accumulation_coef_(0.0),
      sorbed_accumulation_coef_(0.0),
      por_sat_den_vol_(0.0),
      activity_model_(NULL),
      fixed_accumulation_(),
      residual_(),
      prev_molal_(),
      rhs_(),
      jacobian_(),
      lu_solver_() {
  // this is ifdef is breaking the formatting tools
  // #ifdef GLENN
  // , solver(NULL) {
  // #endif
  primary_species_.clear();
  minerals_.clear();
  aqComplexRxns_.clear();
  generalKineticRxns_.clear();
  mineral_rates_.clear();
  ion_exchange_rxns_.clear();
  sorption_isotherm_rxns_.clear();
  total_.clear();
  total_sorbed_.clear();
}  // end Beaker() constructor

Beaker::~Beaker() {
  delete activity_model_;

  if (mineral_rates_.size() != 0) {
    for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
         rate != mineral_rates_.end(); rate++) {
      delete(*rate);
    }
  }

  if (sorption_isotherm_rxns_.size() != 0) {
    std::vector<SorptionIsothermRxn>::iterator rxn;
    for (rxn = sorption_isotherm_rxns_.begin();
         rxn != sorption_isotherm_rxns_.end(); ++rxn) {
      rxn->CleanMemory();
    }
  }

#ifdef GLENN
  delete solver;
#endif
}  // end ~Beaker()

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
}  // end Setup()

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

  parameters.mineral_specific_surface_area.clear();
  parameters.sorption_site_density.clear();
  parameters.isotherm_kd.clear();
  parameters.isotherm_langmuir_b.clear();
  parameters.isotherm_freundlich_n.clear();

  return parameters;
}  // end GetDefaultParameters()

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

  parameters.mineral_specific_surface_area.resize(minerals_.size(), 0.0);
  for (int m = 0; m < minerals_.size(); ++m) {
    double ssa = minerals_.at(m).specific_surface_area();
    parameters.mineral_specific_surface_area.at(m) = ssa;
  }

  if (sorption_isotherm_rxns_.size() > 0) {
    parameters.isotherm_kd.resize(sorption_isotherm_rxns_.size(), 0.0);
    parameters.isotherm_langmuir_b.resize(sorption_isotherm_rxns_.size(), 0.0);
    parameters.isotherm_freundlich_n.resize(sorption_isotherm_rxns_.size(), 0.0);

    std::vector<SorptionIsothermRxn>::const_iterator rxn;
    for (rxn = sorption_isotherm_rxns_.begin(); 
         rxn != sorption_isotherm_rxns_.end(); ++rxn) {
      std::vector<double> params;
      params = rxn->GetIsothermParameters();
      parameters.isotherm_kd.at(rxn->species_id()) = params.at(0);
      if (rxn->IsothermName() == "freundlich") {
        parameters.isotherm_freundlich_n.at(rxn->species_id()) = params.at(1);
      } else if (rxn->IsothermName() == "langmuir") {
        parameters.isotherm_langmuir_b.at(rxn->species_id()) = params.at(1);
      }
    }
  } else {
    parameters.isotherm_kd.clear();
    parameters.isotherm_langmuir_b.clear();
    parameters.isotherm_freundlich_n.clear();
  }

  if (surfaceComplexationRxns_.size() > 0) {
    parameters.sorption_site_density.resize(surfaceComplexationRxns_.size(), 0.0);
    std::vector<SurfaceComplexationRxn>::const_iterator rxn;
    for (rxn = surfaceComplexationRxns_.begin(); 
         rxn != surfaceComplexationRxns_.end(); ++rxn) {
      parameters.sorption_site_density.at(rxn->SiteId()) = rxn->GetSiteDensity();
    }
  } else {
    parameters.sorption_site_density.clear();
  }

  return parameters;
}  // end GetCurrentParameters()

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

  if (parameters.mineral_specific_surface_area.size() > 0) {
    for (int m = 0; m < minerals_.size(); ++m) {
      minerals_.at(m).set_specific_surface_area(parameters.mineral_specific_surface_area.at(m));
    }
  }

  if (parameters.sorption_site_density.size() > 0) {
    for (int scr = 0; scr < surfaceComplexationRxns_.size(); ++scr) {
      surfaceComplexationRxns_.at(scr).UpdateSiteDensity(parameters.sorption_site_density.at(scr));
    }
  }

  if (parameters.isotherm_kd.size() > 0) {
    std::vector<SorptionIsothermRxn>::iterator rxn;
    for (rxn = sorption_isotherm_rxns_.begin();
         rxn != sorption_isotherm_rxns_.end(); ++rxn) {
      std::vector<double> params(2, 0.0);
      int id = rxn->species_id();
      params.at(0) = parameters.isotherm_kd.at(id);
      if (rxn->IsothermName() == "langmuir") {
        params.at(1) = parameters.isotherm_langmuir_b.at(id);
      } else if (rxn->IsothermName() == "freundlich") {
        params.at(1) = parameters.isotherm_freundlich_n.at(id);
      } // else linear --> do nothing
      rxn->SetIsothermParameters(params);
    }
  }

}  // end SetParameters()

//
// public "computation engine" routines
//

// if no water density provided, default is 1000.0 kg/m^3
int Beaker::Speciate(const Beaker::BeakerComponents& components,
                     const Beaker::BeakerParameters& parameters) {
  std::stringstream message;
  double speciation_tolerance = 1.e-12;
  double residual_tolerance = 1.e-12;
  ResetStatus();
  UpdateParameters(parameters, 0.0);
  CheckChargeBalance(components.total);
  // initialize free-ion concentrations, these are actual
  // concentrations, so the value must be positive or ln(free_ion)
  // will return a nan!
  if (components.free_ion.size() > 0) {
    InitializeMolalities(components.free_ion);
  } else {
    InitializeMolalities(1.e-9);
  }

  for (unsigned int m = 0; m < minerals_.size(); m++) {
    minerals_.at(m).set_volume_fraction(components.minerals.at(m));
    // NOTE: SSA already updated by the call to UpdateParameters()!
  }

  // store current molalities
  for (int i = 0; i < ncomp(); i++) {
    prev_molal_.at(i) = primary_species().at(i).molality();
  }

  // TODO(bandre): update CEC from components!

  double max_rel_change;
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
      residual_.at(i) = total_.at(i) - components.total.at(i);
    }
    // add derivatives of total with respect to free to Jacobian
    // units of Jacobian: kg water/sec
    jacobian_.Zero();
    CalculateDTotal();
    jacobian_.AddValues(0, 0, &dtotal_);

    for (int i = 0; i < ncomp(); i++) {
      rhs_.at(i) = residual_.at(i);
    }

    if (debug()) {
      // geh
      message.str("");
      message << "\n- Iteration " << num_iterations << " --------\n";
      chem_out->Write(kDebugBeaker, message);
      DisplayResults();
    }

    if (debug()) {
      print_linear_system("before scale", jacobian_, rhs_);
    }

    // scale the Jacobian
    ScaleRHSAndJacobian();

    if (debug()) {
      print_linear_system("after scale", jacobian_, rhs_);
    }
    // for derivatives with respect to ln concentration, scale columns
    // by primary species concentrations
    for (int i = 0; i < ncomp(); i++) {
      jacobian_.ScaleColumn(i, primary_species().at(i).molality());
    }

    if (debug()) {
      print_linear_system("before solve", jacobian_, rhs_);
    }

    // call solver
    lu_solver_.Solve(&jacobian_, &rhs_);

    if (debug()) {
      print_linear_system("after solve", rhs_);
    }

    // calculate update truncating at a maximum of 5 in log space
    UpdateMolalitiesWithTruncation(5.0);
    // calculate maximum relative change in concentration over all species
    max_rel_change = CalculateMaxRelChangeInMolality();

    if (debug()) {
      message.str("");
      for (int i = 0; i < ncomp(); i++) {
        message << primary_species().at(i).name() << " "
                << primary_species().at(i).molality() << " " << total_.at(i) << "\n";
      }
      chem_out->Write(kDebugBeaker, message);
    }

    num_iterations++;

    max_residual = 0.;
    for (int i = 0; i < ncomp(); i++)
      if (fabs(residual_.at(i)) > max_residual) {
        max_residual = fabs(residual_.at(i));
      }

    // if max_rel_change small enough, turn on activity coefficients
    if (max_rel_change < speciation_tolerance ||
        max_residual < residual_tolerance) {
      calculate_activity_coefs = true;
    }

    // exist if maximum relative change is below tolerance
  } while (max_rel_change > speciation_tolerance &&
           max_residual > residual_tolerance &&
           num_iterations < max_iterations() &&
           !calculate_activity_coefs);

  // for now, initialize total sorbed concentrations based on the current free
  // ion concentrations
  UpdateEquilibriumChemistry();
  status_.num_newton_iterations = num_iterations;
  if (max_rel_change < tolerance()) {
    status_.converged = true;
  }

  if (debug()) {
    message.str("");
    message << "Beaker::speciate: status.num_rhs_evaluations: " << status_.num_rhs_evaluations << std::endl;
    message << "Beaker::speciate: status.num_jacobian_evaluations: " << status_.num_jacobian_evaluations << std::endl;
    message << "Beaker::speciate: status.num_newton_iterations: " << status_.num_newton_iterations << std::endl;
    message << "Beaker::speciate: status.converged: " << status_.converged << std::endl;
    chem_out->Write(kDebugBeaker, message);
  }
  return num_iterations;
}  // end Speciate()

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

  // initialize free-ion concentrations, these are actual
  // concentrations, so the value must be positive or ln(free_ion)
  // will return a nan!
  if (components->free_ion.size() > 0) {
    InitializeMolalities(components->free_ion);
    // testing only    InitializeMolalities(1.e-9);
  }

  // store current molalities
  for (int i = 0; i < ncomp(); i++) {
    prev_molal_.at(i) = primary_species().at(i).molality();
  }

  for (unsigned int m = 0; m < minerals_.size(); m++) {
    minerals_.at(m).set_volume_fraction(components->minerals.at(m));
    // NOTE: SSA already updated by the call to UpdateParameters()!
  }

  // TODO(bandre): update CEC from components!

  // initialize to a large number (not necessary, but safe)
  double max_rel_change = 1.e20;
  // iteration counter
  unsigned int num_iterations = 0;

  // lagging activity coefficients by a time step in this case
  UpdateActivityCoefficients();

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

    if (debug()) {
      print_linear_system("before scale", jacobian_, rhs_);
    }

    // scale the Jacobian
    ScaleRHSAndJacobian();

    if (debug()) {
      print_linear_system("after scale", jacobian_, rhs_);
    }
    // for derivatives with respect to ln concentration, scale columns
    // by primary species concentrations
    for (int i = 0; i < ncomp(); i++) {
      jacobian_.ScaleColumn(i, primary_species().at(i).molality());
    }

    if (debug()) {
      print_linear_system("before solve", jacobian_, rhs_);
    }

    // solve J dlnc = r
    lu_solver_.Solve(&jacobian_, &rhs_);

    if (debug()) {
      print_linear_system("after solve", rhs_);
    }

    // units of solution: mol/kg water (change in molality)
    // calculate update truncating at a maximum of 5 in nat log space
    // update with exp(-dlnc)
    UpdateMolalitiesWithTruncation(5.0);
    // calculate maximum relative change in concentration over all species
    max_rel_change = CalculateMaxRelChangeInMolality();

    if (debug()) {
      message.str("");
      for (int i = 0; i < ncomp(); i++) {
        message << primary_species().at(i).name() << " " <<
            primary_species().at(i).molality() << " " << total_.at(i) << "\n";
      }
      chem_out->Write(kDebugBeaker, message);
    }

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
    error_stream << "Warning: tolerance = " << tolerance() << std::endl;
    error_stream << "Warning: max iterations = " << max_iterations() << std::endl;
    // update before leaving so that we can see the erroneous values!
    UpdateComponents(components);
    Exceptions::amanzi_throw(ChemistryMaxIterationsReached(error_stream.str()));
  }

  status_.num_newton_iterations = num_iterations;
  if (max_rel_change < tolerance()) {
    status_.converged = true;
  }

  // update total concentrations
  UpdateEquilibriumChemistry();
  UpdateKineticMinerals();
  UpdateComponents(components);
  ValidateSolution();

  return num_iterations;
}  // end ReactionStep()


void Beaker::UpdateComponents(Beaker::BeakerComponents* components) {
  // Copy the beaker primary variables into a component struct that
  // can be returned to the driver
  for (int i = 0; i < ncomp(); i++) {
    components->total.at(i) = total_.at(i);
    if (components->free_ion.size() > 0) {
      components->free_ion.at(i) = primary_species().at(i).molality();
    }
    if (components->total_sorbed.size() > 0 &&
        total_sorbed_.size() > 0) {
      components->total_sorbed.at(i) = total_sorbed_.at(i);
    }
  }
  for (unsigned int m = 0; m < minerals_.size(); m++) {
    components->minerals.at(m) = minerals_.at(m).volume_fraction();
  }
  for (unsigned int i = 0; i < ion_exchange_rxns_.size(); i++) {
    components->ion_exchange_sites.at(i) = ion_exchange_rxns_.at(i).site().cation_exchange_capacity();
  }
}  // end UpdateComponents()

void Beaker::CopyComponents(const Beaker::BeakerComponents& from,
                            Beaker::BeakerComponents* to) {
  if (to->total.size() != from.total.size()) {
    to->total.resize(from.total.size());
  }
  for (unsigned int i = 0; i < from.total.size(); i++) {
    to->total.at(i) = from.total.at(i);
  }

  if (to->free_ion.size() != from.free_ion.size()) {
    to->free_ion.resize(from.free_ion.size());
  }
  for (unsigned int i = 0; i < from.free_ion.size(); i++) {
    to->free_ion.at(i) = from.free_ion.at(i);
  }

  if (to->total_sorbed.size() != from.total_sorbed.size()) {
    to->total_sorbed.resize(from.total_sorbed.size());
  }
  for (unsigned int i = 0; i < from.total_sorbed.size(); i++) {
    to->total_sorbed.at(i) = from.total_sorbed.at(i);
  }

  if (to->minerals.size() != from.minerals.size()) {
    to->minerals.resize(from.minerals.size());
  }
  for (unsigned int i = 0; i < from.minerals.size(); i++) {
    to->minerals.at(i) = from.minerals.at(i);
  }

  if (to->ion_exchange_sites.size() != from.ion_exchange_sites.size()) {
    to->ion_exchange_sites.resize(from.ion_exchange_sites.size());
  }
  for (unsigned int i = 0; i < from.ion_exchange_sites.size(); i++) {
    to->ion_exchange_sites.at(i) = from.ion_exchange_sites.at(i);
  }
}  // end Beaker::CopyComponents()

void Beaker::GetPrimaryNames(std::vector<std::string>* names) const {
  names->clear();
  for (std::vector<Species>::const_iterator primary = primary_species().begin();
       primary != primary_species().end(); primary++) {
    names->push_back(primary->name());
  }
}  // end GetPrimaryNames()

int Beaker::GetPrimaryIndex(const std::string& name) const {
  int index = -1;
  for (std::vector<Species>::const_iterator primary = primary_species().begin();
       primary != primary_species().end(); primary++) {
    if (primary->name() == name) {
      index = primary->identifier();
    }
  }
  return index;
}  // end GetPrimaryIndex()

//
// public output
//

void Beaker::Display(void) const {
  chem_out->Write(kVerbose, "-- Beaker description ------------------------------------------------\n");
  DisplayParameters();

  DisplayPrimary();

  DisplayAqueousEquilibriumComplexes();

  DisplayGeneralKinetics();

  DisplayMinerals();

  DisplayMineralKinetics();

  DisplayIonExchangeSites();

  DisplayIonExchangeComplexes();

  DisplaySurfaceSites();

  DisplaySurfaceComplexes();

  DisplaySorptionIsotherms();

  chem_out->Write(kVerbose, "------------------------------------------------ Beaker description --\n");
}  // end Display()

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
              << components.minerals.at(m) << std::endl;
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
  chem_out->Write(kVerbose, message);
}  // end DisplayComponents

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
  chem_out->Write(kVerbose, message);

  primary_species().at(0).DisplayResultsHeader();
  for (int i = 0; i < ncomp(); i++) {
    primary_species().at(i).DisplayResults();
  }

  // same header info as primaries....
  for (unsigned int i = 0; i < aqComplexRxns_.size(); i++) {
    aqComplexRxns_.at(i).DisplayResults();
  }

  if (minerals_.size() > 0) {
    chem_out->Write(kVerbose, "---- Minerals\n");
    minerals_[0].DisplayResultsHeader();
    for (unsigned int i = 0; i < minerals_.size(); i++) {
      minerals_.at(i).DisplayResults();
    }
  }

  if (ion_exchange_rxns_.size() > 0) {
    chem_out->Write(kVerbose, "---- Ion Exchange Sites\n");
    ion_exchange_rxns_[0].DisplayResultsHeader();
    for (unsigned int i = 0; i < ion_exchange_rxns_.size(); i++) {
      ion_exchange_rxns_.at(i).site().DisplayResults();
    }
  }

  if (ion_exchange_rxns_.size() > 0) {
    chem_out->Write(kVerbose, "---- Ion Exchange Complexes\n");
    ion_exchange_rxns_[0].DisplayResultsHeader();
    for (unsigned int i = 0; i < ion_exchange_rxns_.size(); i++) {
      for (unsigned int j = 0; j < ion_exchange_rxns_.at(i).ionx_complexes().size(); j++) {
        (ion_exchange_rxns_.at(i).ionx_complexes())[j].DisplayResults();
      }
    }
  }

  if (surfaceComplexationRxns_.size() > 0) {
    chem_out->Write(kVerbose, "---- Surface Complexation Reactions\n");
    for (unsigned int i = 0; i < surfaceComplexationRxns_.size(); i++) {
      surfaceComplexationRxns_.at(i).DisplayResultsHeader();
      surfaceComplexationRxns_.at(i).DisplayResults();
    }
  }

  chem_out->Write(kVerbose, "---------------------------------------------------------- Solution --\n\n");
}  // end DisplayResults()

void Beaker::DisplayTotalColumnHeaders(void) const {
  std::stringstream message;
  message << std::setw(15) << "Time (s)";
  for (int i = 0; i < ncomp(); i++) {
    message << std::setw(15) << primary_species().at(i).name();
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
  chem_out->Write(kVerbose, message);
}  // end DisplayTotalColumnHeaders()

void Beaker::DisplayTotalColumns(const double time,
                                 const BeakerComponents& components) const {
  std::stringstream message;
  message << std::scientific << std::setprecision(6) << std::setw(15);
  message << time;
  for (int i = 0; i < ncomp(); i++) {
    message << std::setw(15) << components.total.at(i);
  }
  if (total_sorbed_.size() > 0) {
    for (int i = 0; i < total_sorbed_.size(); i++) {
      message << std::setw(15) << components.total_sorbed.at(i);
    }
  }
  if (minerals().size() > 0) {
    for (int m = 0; m < minerals().size(); ++m) {
      message << std::setw(15) << components.minerals.at(m);
    }
  }
  message << std::endl;
  chem_out->Write(kVerbose, message);
}  // end DisplayTotalColumns()




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
      sorption_isotherm_rxns_.size() > 0) {
    total_sorbed_.resize(ncomp(), 0.0);
    dtotal_sorbed_.Resize(ncomp());
    dtotal_sorbed_.Zero();
  } else {
    total_sorbed_.resize(0);
    //dtotal_sorbed_.Resize(0);
  }

  fixed_accumulation_.resize(ncomp());
  residual_.resize(ncomp());
  prev_molal_.resize(ncomp());

  jacobian_.Resize(ncomp());
  rhs_.resize(ncomp());
  lu_solver_.Initialize(ncomp());
}  // end ResizeInternalMemry()



void Beaker::VerifyComponentSizes(const Beaker::BeakerComponents& components) const {
  // some helpful error checking goes here...
  bool error = false;
  std::ostringstream error_stream;
  error_stream << "Beaker::VerifyComponentSizes(): database input and component initial conditions do not match:\n";

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

  if (this->minerals().size() != components.minerals.size()) {
    error = true;
    error_stream << "minerals.size(" << this->minerals().size()
                 << ") and components.minerals.size(" << components.minerals.size()
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
    error_stream << "total_sorbed.size and components.total_sorbed.size do not match.\n";
  }

  if (error) {
    Exceptions::amanzi_throw(ChemistryMemorySizeError(error_stream.str()));
  }
}  // end VerifyComponentSizes()

void Beaker::SetComponents(const Beaker::BeakerComponents& components) {
  // TODO(bandre): Not sure this is actually doing anything as currently
  // used. Should be renamed ApplyComponentsToBeaker() or something
  // similar, then repurposed to handle this behavior uniformly in the
  // SimpleThermoDatabase::Setup(), ReactionStep() an Speciate() calls.
  unsigned int size = components.ion_exchange_sites.size();
  if (ion_exchange_rxns().size() == size) {
    // for now, we will set the size of total_sorbed = ncomp() if # sites > 0
    for (unsigned int ies = 0; ies < size; ies++) {
      ion_exchange_rxns_[ies].set_cation_exchange_capacity(components.ion_exchange_sites.at(ies));
    }
  } else {
    std::ostringstream error_stream;
    error_stream << "Beaker::SetComponents(): \n";
    error_stream << "ion_exchange_sites.size and components.ion_exchange_sites.size do not match.\n";
    Exceptions::amanzi_throw(ChemistryUnrecoverableError(error_stream.str()));
  }

  size = components.minerals.size();
  if (minerals().size() == size) {
    for (unsigned int m = 0; m < size; m++) {
      minerals_.at(m).set_volume_fraction(components.minerals.at(m));
    }
  } else {
    std::ostringstream error_stream;
    error_stream << "Beaker::SetComponents(): \n";
    error_stream << "minerals.size and components.minerals.size do not match.\n";
    Exceptions::amanzi_throw(ChemistryUnrecoverableError(error_stream.str()));
  }
}  // end SetComponents()

void Beaker::SetupActivityModel(std::string model,
                                std::string pitzer_database,
                                std::string jfunction_pitzer) {

  delete activity_model_;

  ActivityModel::ActivityModelParameters parameters;
  parameters.database_filename = pitzer_database;
  parameters.pitzer_jfunction = jfunction_pitzer;
  if (verbosity() == kDebugActivityModel) {
    parameters.verbosity = verbosity();
  }

  ActivityModelFactory amf;

  activity_model_ = amf.Create(model, parameters,
                               primary_species(), aqComplexRxns_);

  if (debug()) {
    activity_model_->Display();
  }
}  // end SetupActivityModel()

void Beaker::AddPrimarySpecies(const Species& s) {
  primary_species_.push_back(s);
}  // end AddPrimarySpecies()

void Beaker::AddIonExchangeRxn(const IonExchangeRxn& ionx_rxn) {
  ion_exchange_rxns_.push_back(ionx_rxn);
}  // end AddIonExchangeRxn()

void Beaker::AddIonExchangeComplex(const int irxn, const IonExchangeComplex& ionx_complex) {
  ion_exchange_rxns_[irxn].AddIonExchangeComplex(ionx_complex);
}  // end AddIonExchangeRxn()

void Beaker::AddAqueousEquilibriumComplex(const AqueousEquilibriumComplex& c) {
  aqComplexRxns_.push_back(c);
}  // end AddAqueousEquilibriumComplex()

void Beaker::AddMineral(const Mineral& m) {
  minerals_.push_back(m);
}  // end AddMineral()

void Beaker::AddMineralKineticRate(KineticRate* rate) {
  mineral_rates_.push_back(rate);
}  // end AddMineralKineticRate()

bool Beaker::HaveKinetics(void) const {
  bool have_kinetics = false;
  if (mineral_rates_.size()) {
    have_kinetics = true;
  }
  // add other kinetic processes here....

  return have_kinetics;
}  // end HaveKinetics()

void Beaker::AddGeneralRxn(const GeneralRxn& r) {
  generalKineticRxns_.push_back(r);
}  // end AddGeneralRxn()

void Beaker::AddSurfaceComplexationRxn(const SurfaceComplexationRxn& r) {
  surfaceComplexationRxns_.push_back(r);
}  // end AddSurfaceComplexationRxn()

void Beaker::AddSorptionIsothermRxn(const SorptionIsothermRxn& r) {
  sorption_isotherm_rxns_.push_back(r);
}  // end AddSorptionIsothermRxn()



/*******************************************************************************
 **
 ** private routines
 **
 ******************************************************************************/
void Beaker::UpdateParameters(const Beaker::BeakerParameters& parameters,
                              double delta_t) {
  dt(delta_t);  // delta time = [sec]
  SetParameters(parameters);
}  // end UpdateParameters()

void Beaker::ResetStatus(void) {
  status_.num_rhs_evaluations = 0;
  status_.num_jacobian_evaluations = 0;
  status_.num_newton_iterations = 0;
  status_.converged = false;
}

void Beaker::update_accumulation_coefficients(void) {
  aqueous_accumulation_coef(porosity() * saturation() * volume() * 1000.0 / dt());
  sorbed_accumulation_coef(volume() / dt());
}  // end update_accumulation_coefficients

void Beaker::update_por_sat_den_vol(void) {
  por_sat_den_vol(porosity() * saturation() * water_density_kg_m3() * volume());
}  // end update_por_sat_den_vol()


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
}  // end UpdateActivityCoefficients

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

}  // end UpdateKineticMinerals()

void Beaker::InitializeMolalities(double initial_molality) {
  for (std::vector<Species>::iterator i = primary_species_.begin();
       i != primary_species_.end(); i++) {
    i->update(initial_molality);
  }
}  // end InitializeMolalities

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
}  // end InitializeMolalities()

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
}  // end UpdateEquilibriumChemistry()

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

}  // end CalculateTotal()

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
}  // end CalculateDTotal()

void Beaker::UpdateKineticChemistry(void) {
  // loop over general kinetic reactions and update effective rates
  for (std::vector<GeneralRxn>::iterator i = generalKineticRxns_.begin();
       i != generalKineticRxns_.end(); i++) {
    i->update_rates(primary_species());
  }

  // add mineral saturation and rate calculations here
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->Update(primary_species(), minerals_);
  }
  // add multirate kinetic surface complexation reaction quotient calculations
  // here
}  // end UpdateKineticChemistry()

void Beaker::AddKineticChemistryToResidual(void) {
  // loop over general kinetic reactions and add rates
  for (std::vector<GeneralRxn>::iterator i = generalKineticRxns_.begin();
       i != generalKineticRxns_.end(); i++) {
    i->addContributionToResidual(&residual_, por_sat_den_vol());
  }

  // add mineral mineral contribution to residual here.  units = mol/sec.
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->AddContributionToResidual(minerals_, volume(), &residual_);
  }

  // add multirate kinetic surface complexation contribution to residual here.
}  // end AddKineticChemistryToResidual()

void Beaker::AddKineticChemistryToJacobian(void) {
  // loop over general kinetic reactions and add rates
  for (std::vector<GeneralRxn>::iterator i = generalKineticRxns_.begin();
       i != generalKineticRxns_.end(); i++) {
    i->addContributionToJacobian(&jacobian_, primary_species(), por_sat_den_vol());
  }

  // add mineral mineral contribution to Jacobian here.  units = kg water/sec.
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->AddContributionToJacobian(primary_species(), minerals_, volume(), &jacobian_);
  }

  // add multirate kinetic surface complexation contribution to Jacobian here.
}  // end AddKineticChemistryToJacobian()



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
}  // end AddAccumulation()

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
}  // end AddAccumulationDerivative()

void Beaker::CalculateFixedAccumulation(const std::vector<double>& total,
                                        const std::vector<double>& total_sorbed,
                                        std::vector<double>* fixed_accumulation) {
  for (unsigned int i = 0; i < total.size(); i++) {
    fixed_accumulation->at(i) = 0.0;
  }
  AddAccumulation(total, total_sorbed, fixed_accumulation);
}  // end CalculateFixedAccumulation()

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
}  // end CalculateResidual()

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
}  // end CalculateJacobian()

void Beaker::ScaleRHSAndJacobian(void) {
  for (int i = 0; i < jacobian_.size(); i++) {
    double max = jacobian_.GetRowAbsMax(i);
    if (max > 1.0) {
      double scale = 1.0 / max;
      rhs_.at(i) *= scale;
      jacobian_.ScaleRow(i, scale);
    }
  }
}  // end ScaleRHSAndJacobian()

void Beaker::UpdateMolalitiesWithTruncation(double max_change) {
  // truncate the rhs to max_change
  for (int i = 0; i < ncomp(); i++) {
    if (rhs_.at(i) > max_change) {
      rhs_.at(i) = max_change;
    } else if (rhs_.at(i) < -max_change) {
      rhs_.at(i) = -max_change;
    }
    // store the previous solution
    prev_molal_.at(i) = primary_species().at(i).molality();
    // update primary species molalities (log formulation)
    double molality = prev_molal_.at(i) * std::exp(-rhs_.at(i));
    primary_species_.at(i).update(molality);
  }
}  // end UpdateMolalitiesWithTruncation()

double Beaker::CalculateMaxRelChangeInMolality(void) {
  double max_rel_change = 0.0;
  for (int i = 0; i < ncomp(); i++) {
    double delta = fabs(primary_species().at(i).molality() - prev_molal().at(i)) / 
        prev_molal().at(i);
    max_rel_change = delta > max_rel_change ? delta : max_rel_change;
  }
  return max_rel_change;
}  // end CalculateMaxRelChangeInMolality()


void Beaker::CheckChargeBalance(const std::vector<double>& aqueous_totals) const {
  double charge_balance = 0.0;
  for (unsigned int i = 0; i < aqueous_totals.size(); i++) {
    charge_balance += aqueous_totals.at(i) * primary_species().at(i).charge();
  }
  if (std::fabs(charge_balance) > tolerance()) {
    std::stringstream message;
    message << "WARNING: Beaker::CheckChargeBalance() :\n"
            << "         charge balance = " << std::scientific
            << charge_balance << std::fixed << std::endl;
    chem_out->Write(kWarning, message);
  }
}  // end CheckChargeBalance()


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
}  // end ValidateSolution()




/*******************************************************************************
 **
 **  Output related functions
 **
 *******************************************************************************/
void Beaker::display(void) const {
  chem_out->Write(kVerbose, "----- Beaker description ------\n");
  chem_out->Write(kVerbose, "Primary Species:");
  for (std::vector<Species>::const_iterator primary = primary_species().begin();
       primary != primary_species().end(); primary++) {
    primary->display();
  }
  chem_out->Write(kVerbose, "\n");
  chem_out->Write(kVerbose, "Aqueous Equilibrium Complexes:\n");
  for (std::vector<AqueousEquilibriumComplex>::const_iterator aec = aqComplexRxns_.begin();
       aec != aqComplexRxns_.end(); aec++) {
    aec->display();
  }
  chem_out->Write(kVerbose, "-------------------------------------");
}  // end display()

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
  chem_out->Write(kVerbose, message);
}  // end DisplayParameters()

void Beaker::DisplayPrimary(void) const {
  std::stringstream message;
  message << "---- Primary Species" << std::endl;
  message << std::setw(15) << "Species"
          << std::setw(10) << "Charge"
          << std::setw(10) << "GMW"
          << std::setw(10) << "D-H a0"
          << std::endl;
  chem_out->Write(kVerbose, message);
  for (std::vector<Species>::const_iterator primary = primary_species().begin();
       primary != primary_species().end(); primary++) {
    primary->Display();
  }
  chem_out->Write(kVerbose, "\n");
}  // end DisplayPrimary()

void Beaker::DisplayAqueousEquilibriumComplexes(void) const {
  std::stringstream message;
  message << "---- Aqueous Equilibrium Complexes" << std::endl;
  message << std::setw(12) << "Reaction"
          << std::setw(38) << "log Keq"
          << std::setw(8) << "Charge"
          << std::setw(10) << "GMW"
          << std::setw(8) << "D-H a0"
          << std::endl;
  chem_out->Write(kVerbose, message);
  for (std::vector<AqueousEquilibriumComplex>::const_iterator aec = aqComplexRxns_.begin();
       aec != aqComplexRxns_.end(); aec++) {
    aec->Display();
  }
  chem_out->Write(kVerbose, "\n");
}  // end DisplayAqueousEquilibriumComplexes()

void Beaker::DisplayGeneralKinetics(void) const {
  if (generalKineticRxns_.size() > 0) {
    std::stringstream message;
    message << "---- General Kinetics" << std::endl;
    message << std::setw(12) << "Reaction" << std::endl;
    chem_out->Write(kVerbose, message);
    for (std::vector<GeneralRxn>::const_iterator rxn = generalKineticRxns_.begin();
         rxn != generalKineticRxns_.end(); rxn++) {
      rxn->Display();
    }
    chem_out->Write(kVerbose, "\n");
  }
}  // end DisplayAqueousEquilibriumComplexes()

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
    chem_out->Write(kVerbose, message);
    for (std::vector<Mineral>::const_iterator m = minerals_.begin();
         m != minerals_.end(); m++) {
      m->Display();
    }
    chem_out->Write(kVerbose, "\n");
  }
}  // end DisplayMinerals()

void Beaker::DisplayMineralKinetics(void) const {
  if (mineral_rates_.size() > 0) {
    std::stringstream message;
    chem_out->Write(kVerbose, "---- Mineral Kinetics\n");
    for (std::vector<KineticRate*>::const_iterator m = mineral_rates_.begin();
         m != mineral_rates_.end(); m++) {
      (*m)->Display();
    }
    chem_out->Write(kVerbose, "\n");
  }
}  // end DisplayMineralKinetics()

void Beaker::DisplayIonExchangeSites(void) const {
  if (ion_exchange_rxns_.size() > 0) {
    std::stringstream message;
    message << "---- Ion Exchange Sites" << std::endl;
    message << std::setw(15) << "Species"
            << std::setw(20) << "Location"
            << std::setw(10) << "Charge"
            << std::setw(10) << "CEC"
            << std::endl;
    chem_out->Write(kVerbose, message);
    std::vector<IonExchangeRxn>::const_iterator rxn;
    for (rxn = ion_exchange_rxns_.begin();
         rxn != ion_exchange_rxns_.end(); rxn++) {
      rxn->site().Display();
    }
    chem_out->Write(kVerbose, "\n");
  }
}  // end DisplayIonExchangeSites()

void Beaker::DisplayIonExchangeComplexes(void) const {
  if (ion_exchange_rxns_.size() > 0) {
    std::stringstream message;
    message << "---- Ion Exchange Complexes" << std::endl;
    message << std::setw(12) << "Reaction"
            << std::setw(38) << "K"
            << std::endl;
    chem_out->Write(kVerbose, message);
    std::vector<IonExchangeRxn>::const_iterator ier;
    for (ier = ion_exchange_rxns_.begin();
         ier != ion_exchange_rxns_.end(); ier++) {
      ier->Display();
    }
    chem_out->Write(kVerbose, "\n");
  }
}  // end DisplayIonExchangeComplexes()

void Beaker::DisplaySurfaceSites(void) const {
  if (total_sorbed_.size() > 0) {
    std::stringstream message;
    message << "---- Surface Sites" << std::endl;
    message << std::setw(15) << "Species"
            << std::setw(15) << "Site Density"
            << std::endl;
    chem_out->Write(kVerbose, message);
    std::vector<SurfaceComplexationRxn>::const_iterator s;
    for (s = surfaceComplexationRxns_.begin();
         s != surfaceComplexationRxns_.end(); s++) {
      s->DisplaySite();
    }
    chem_out->Write(kVerbose, "\n");
  }
}  // end DisplaySurfaceSites()

void Beaker::DisplaySurfaceComplexes(void) const {
  if (total_sorbed_.size() > 0) {
    std::stringstream message;
    message << "---- Surface Complexes" << std::endl;
    message << std::setw(12) << "Reaction"
            << std::setw(38) << "log Keq"
            << std::setw(10) << "charge"
            << std::endl;
    chem_out->Write(kVerbose, message);
    std::vector<SurfaceComplexationRxn>::const_iterator s;
    for (s = surfaceComplexationRxns_.begin();
         s != surfaceComplexationRxns_.end(); s++) {
      s->DisplayComplexes();
    }
    chem_out->Write(kVerbose, "\n");
  }
}  // end DisplaySurfaceComplexes()

void Beaker::DisplaySorptionIsotherms(void) const {
  if (total_sorbed_.size() > 0) {
    std::stringstream message;
    message << "---- Equilibrium Sorption Isotherms" << std::endl;
    message << std::setw(12) << "Species"
            << std::setw(15) << "isotherm"
            << std::setw(15) << "parameters"
            << std::endl;
    chem_out->Write(kVerbose, message);
    std::vector<SorptionIsothermRxn>::const_iterator s;
    for (s = sorption_isotherm_rxns_.begin();
         s != sorption_isotherm_rxns_.end(); s++) {
      s->Display();
    }
    chem_out->Write(kVerbose, "\n");
  }
}  // end DisplaySurfaceComplexes()


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
}  // end print_results()

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
}  // end print_results()

void Beaker::print_linear_system(const std::string& s, 
                                 const MatrixBlock& A,
                                 const std::vector<double>& vector) const {
  std::stringstream message;
  message << s << std::endl;
  for (unsigned int i = 0; i < vector.size(); i++) {
    message << "RHS: " << primary_species().at(i).name() << " " << vector.at(i) << std::endl;
  }
  chem_out->Write(kVerbose,message);
  A.Print_ij();
}  // end print_linear_system()

void Beaker::print_linear_system(const std::string& s, 
                                 const std::vector<double>& vector) const {
  std::stringstream message;
  message << s << std::endl;
  for (unsigned int i = 0; i < vector.size(); i++) {
    message << "RHS: " << primary_species().at(i).name() << " " << vector.at(i) << std::endl;
  }
  chem_out->Write(kVerbose,message);
}  // end print_linear_system()

}  // namespace chemistry
}  // namespace amanzi
