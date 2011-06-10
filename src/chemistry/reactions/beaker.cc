/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** TODO: update mineral volume fractions to components.minerals after kinetics....
**
** TODO: finish implementing ion exchange jacobian
*/

#include "beaker.hh"

#include <cstdlib>

#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>

#include "activity_model.hh"
#include "activity_model_factory.hh"
#include "aqueous_equilibrium_complex.hh"
#include "block.hh"
#include "general_rxn.hh"
#include "chemistry_exception.hh"
#include "ion_exchange_complex.hh"
#include "ion_exchange_site.hh"
#include "kinetic_rate.hh"
#include "lu.hh"
#include "mineral.hh"
#include "mineral_kinetics_factory.hh"
#include "species.hh"
#include "surface_complexation_rxn.hh"
#include "verbosity.hh"

#include "exceptions.hh"

// solver defaults
const double Beaker::tolerance_default = 1.0e-12;
const unsigned int Beaker::max_iterations_default = 250;
// default physical parameters
const double Beaker::porosity_default = 1.0;  // [-]
const double Beaker::saturation_default = 1.0;  // [-]
const double Beaker::water_density_kg_m3_default = 1000.0;
const double Beaker::volume_default = 1.0;  // [m^3]

Beaker::Beaker()
    : verbosity_(kSilent),
      tolerance_(tolerance_default),
      max_iterations_(max_iterations_default),
      ncomp_(0),
      dtotal_(NULL),
      dtotal_sorbed_(NULL),
      porosity_(porosity_default),
      saturation_(saturation_default),
      water_density_kg_m3_(water_density_kg_m3_default),
      water_density_kg_L_(1.0),
      volume_(volume_default),
      dt_(0.0),
      aqueous_accumulation_coef_(0.0),
      sorbed_accumulation_coef_(0.0),
      por_sat_den_vol_(0.0),
      activity_model_(NULL),
      J(NULL)
#ifdef GLENN
    , solver(NULL)
#endif
{
  primarySpecies_.clear();
  minerals_.clear();
  ion_exchange_sites_.clear();
  aqComplexRxns_.clear();
  generalKineticRxns_.clear();
  mineral_rates_.clear();
  ion_exchange_rxns_.clear();

  total_.clear();
  total_sorbed_.clear();
}  // end Beaker() constructor

Beaker::~Beaker() {
  // According to C++ In A Nutshell (and other sources), no check
  // on NULL is necessary for "delete".
  delete dtotal_;
  delete dtotal_sorbed_;
  delete J;
  delete activity_model_;

  if (mineral_rates_.size() != 0) {
    for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
         rate != mineral_rates_.end(); rate++) {
      delete(*rate);
    }
  }
#ifdef GLENN
  delete solver;
#endif
}  // end Beaker destructor

void Beaker::resize() {
  delete dtotal_;
  dtotal_ = new Block(ncomp());
  // For now, assume that dtotal is always allocated
  delete dtotal_sorbed_;
  dtotal_sorbed_ = new Block(ncomp());
  dtotal_sorbed_->zero();

  fixed_accumulation.resize(ncomp());
  residual.resize(ncomp());
  prev_molal.resize(ncomp());
  total_.resize(ncomp());
  // TODO: FIXED? this should only be done if we are actually using sorption.
  if (surfaceComplexationRxns_.size() > 0) {
    total_sorbed_.resize(ncomp(), 0.0);
    // for (unsigned int i = 0; i < total_sorbed_.size(); i++)
    //  total_sorbed_[i] = 0.;
  }

  delete J;
  J = new Block(ncomp());
  rhs.resize(ncomp());
  indices.resize(ncomp());

#ifdef GLENN
  solver = new DirectSolver();
  solver->Initialize(ncomp());
#endif
}  // end resize()

void Beaker::resize(int n) {
  ncomp(n);
  resize();
}  // end resize()

void Beaker::Setup(const Beaker::BeakerComponents& components,
                   const Beaker::BeakerParameters& parameters) {
  SetParameters(parameters);

  this->SetupActivityModel(parameters.activity_model_name);
  this->resize((int)this->primary_species().size());
  this->VerifyComponentSizes(components);
}  // end Setup()

void Beaker::VerifyComponentSizes(const Beaker::BeakerComponents& components) {
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

  if (this->ion_exchange_sites().size() != components.ion_exchange_sites.size()) {
    error = true;
    error_stream << "ion_exchange_sites.size(" << this->ion_exchange_sites().size()
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
    Exceptions::amanzi_throw(ChemistryUnrecoverableError(error_stream.str()));
  }

}  // end VerifyComponentSizes()

void Beaker::SetComponents(const Beaker::BeakerComponents& components) {
  unsigned int size = components.ion_exchange_sites.size();
  if (ion_exchange_sites().size() == size) {
    for (unsigned int ies = 0; ies < size; ies++) {
      ion_exchange_sites_[ies].set_cation_exchange_capacity(components.ion_exchange_sites.at(ies));
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
      minerals_.at(m).UpdateSurfaceAreaFromVolumeFraction(volume());
    }
  } else {
    std::ostringstream error_stream;
    error_stream << "Beaker::SetComponents(): \n";
    error_stream << "minerals.size and components.minerals.size do not match.\n";
    Exceptions::amanzi_throw(ChemistryUnrecoverableError(error_stream.str()));
  }
}  // end SetComponents()

void Beaker::SetupActivityModel(std::string model) {
  if (model != ActivityModelFactory::unit &&
      model != ActivityModelFactory::debye_huckel) {
    model = ActivityModelFactory::unit;
  }
  if (activity_model_ != NULL) {
    delete activity_model_;
  }
  ActivityModelFactory amf;

  activity_model_ = amf.Create(model);
  if (verbosity() == kDebugBeaker) {
    activity_model_->Display();
  }
}  // end SetupActivityModel()

void Beaker::addPrimarySpecies(const Species& s) {
  primarySpecies_.push_back(s);
}  // end addPrimarySpecies()

void Beaker::AddIonExchangeSite(const IonExchangeSite& exchanger) {
  ion_exchange_sites_.push_back(exchanger);
}  // end AddIonExchangeSites()

void Beaker::AddIonExchangeComplex(const IonExchangeComplex& exchange_complex) {
  ion_exchange_rxns_.push_back(exchange_complex);
}  // end AddIonExchangeSites()

void Beaker::addAqueousEquilibriumComplex(const AqueousEquilibriumComplex& c) {
  aqComplexRxns_.push_back(c);
}  // end addAqueousEquilibriumComplex()

void Beaker::addMineral(const Mineral& m) {
  minerals_.push_back(m);
}  // end addMineral()

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

void Beaker::addGeneralRxn(const GeneralRxn& r) {
  generalKineticRxns_.push_back(r);
}  // end addGeneralRxn()

void Beaker::addSurfaceComplexationRxn(const SurfaceComplexationRxn& r) {
  surfaceComplexationRxns_.push_back(r);
}  // end addSurfaceComplexationRxn()

Beaker::BeakerParameters Beaker::GetDefaultParameters(void) const {
  Beaker::BeakerParameters parameters;

  parameters.thermo_database_file.clear();

  parameters.tolerance = tolerance_default;
  parameters.max_iterations = max_iterations_default;

  parameters.activity_model_name = ActivityModelFactory::unit;

  parameters.porosity = porosity_default;
  parameters.saturation = saturation_default;
  parameters.water_density = water_density_kg_m3_default;  // kg / m^3
  parameters.volume = volume_default;  // m^3

  return parameters;
}  // end GetDefaultParameters()

Beaker::BeakerParameters Beaker::GetCurrentParameters(void) const {
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
}  // end GetCurrentParameters()

void Beaker::SetParameters(const Beaker::BeakerParameters& parameters) {
  tolerance(parameters.tolerance);
  max_iterations(parameters.max_iterations);
  porosity(parameters.porosity);
  water_density_kg_m3(parameters.water_density);  // den = [kg/m^3]
  saturation(parameters.saturation);
  volume(parameters.volume);  // vol = [m^3]
  update_accumulation_coefficients();
  update_por_sat_den_vol();
}  // end SetParameters()


void Beaker::updateParameters(const Beaker::BeakerParameters& parameters,
                              double delta_t) {
  dt(delta_t);  // delta time = [sec]
  SetParameters(parameters);
}  // end updateParameters()

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


void Beaker::updateActivityCoefficients() {

  // return;
  activity_model_->CalculateIonicStrength(primarySpecies_,
                                          aqComplexRxns_);
  activity_model_->CalculateActivityCoefficients(&primarySpecies_,
                                                 &aqComplexRxns_);
  for (std::vector<Species>::iterator i = primarySpecies_.begin();
       i != primarySpecies_.end(); i++) {
    i->update();
  }

}

void Beaker::initializeMolalities(double initial_molality) {
  for (std::vector<Species>::iterator i = primarySpecies_.begin();
       i != primarySpecies_.end(); i++) {
    i->update(initial_molality);
  }
}  // end initializeMolalities

void Beaker::initializeMolalities(const std::vector<double>& initial_molalities) {
  if (initial_molalities.size() != primarySpecies_.size()) {
    std::ostringstream error_stream;
    error_stream << "Beaker::initializeMolalities(): \n";
    error_stream << "   Mismatch in size of initial_molalities array "
                 << "and number of primarySpecies" << std::endl;
    Exceptions::amanzi_throw(ChemistryUnrecoverableError(error_stream.str()));
  }

  // iterator doesnt seem to work then passing a vector entry - geh
  for (int i = 0; i < (int)primarySpecies_.size(); i++) {
    primarySpecies_[i].update(initial_molalities[i]);
  }
}  // end initializeMolalities()

void Beaker::updateEquilibriumChemistry(void) {

  //    calculateActivityCoefficients(-1);

  // update primary species activities
  for (std::vector<Species>::iterator primary = primarySpecies_.begin();
       primary != primarySpecies_.end(); primary++) {
    primary->update();
  }
  // calculated seconday aqueous complex concentrations
  for (std::vector<AqueousEquilibriumComplex>::iterator aqcplx =
           aqComplexRxns_.begin();
       aqcplx != aqComplexRxns_.end(); aqcplx++) {
    aqcplx->Update(primarySpecies_);
  }

  // calculate mineral saturation states
  for (std::vector<Mineral>::iterator m = minerals_.begin();
       m != minerals_.end(); m++) {
    m->Update(primarySpecies_);
  }

  // add equilibrium surface complexation here
  for (std::vector<SurfaceComplexationRxn>::iterator srfcplx =
           surfaceComplexationRxns_.begin();
       srfcplx != surfaceComplexationRxns_.end(); srfcplx++) {
    srfcplx->Update(primarySpecies_);
  }

  // add equilibrium ion exchange here?
  for (std::vector<IonExchangeComplex>::iterator iec = ion_exchange_rxns_.begin();
       iec != ion_exchange_rxns_.end(); iec++) {
    iec->Update(primary_species(), ion_exchange_sites());
  }

  // calculate total component concentrations
  calculateTotal();

}  // end updateEquilibriumChemistry()

void Beaker::calculateTotal(std::vector<double> *total,
                            std::vector<double> *total_sorbed) {
  // add in primaries
  for (unsigned int i = 0; i < total->size(); i++) {
    (*total)[i] = primarySpecies_.at(i).molality();
  }

  // add in aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::iterator i = aqComplexRxns_.begin();
       i != aqComplexRxns_.end(); i++) {
    i->AddContributionToTotal(total);
  }

  // scale by water density to convert to molarity
  for (unsigned int i = 0; i < total->size(); i++) {
    (*total)[i] *= water_density_kg_L();
  }
  // calculate sorbed totals
  // initialize to zero
  for (unsigned int i = 0; i < total_sorbed->size(); i++) {
    (*total_sorbed)[i] = 0.;
  }
  // add in contributions
  for (std::vector<SurfaceComplexationRxn>::iterator i =
           surfaceComplexationRxns_.begin();
       i != surfaceComplexationRxns_.end(); i++) {
    i->AddContributionToTotal(total_sorbed);
  }

}  // end calculateTotal()

void Beaker::calculateTotal(void) {
  calculateTotal(&total_, &total_sorbed_);
}  // end calculateTotal()

void Beaker::calculateDTotal(Block* dtotal, Block* dtotal_sorbed) {
  dtotal->zero();
  // derivative with respect to free-ion is 1.
  dtotal->setDiagonal(1.);

  // add in derviative of complex contribution with respect to free-ion
  for (std::vector<AqueousEquilibriumComplex>::iterator i = aqComplexRxns_.begin();
       i != aqComplexRxns_.end(); i++) {
    i->AddContributionToDTotal(primarySpecies_, dtotal);
  }

  // scale by density of water
  dtotal->scale(water_density_kg_L());

  if (dtotal_sorbed) {
    dtotal_sorbed->zero();
  }
  // calculate sorbed derivatives
  for (std::vector<SurfaceComplexationRxn>::iterator i =
           surfaceComplexationRxns_.begin();
       i != surfaceComplexationRxns_.end(); i++) {
    i->AddContributionToDTotal(primarySpecies_, dtotal_sorbed);
  }

}  // end calculateDTotal()

void Beaker::calculateDTotal(void) {
  calculateDTotal(dtotal_, dtotal_sorbed_);
}

void Beaker::updateKineticChemistry(void) {
  // loop over general kinetic reactions and update effective rates
  for (std::vector<GeneralRxn>::iterator i = generalKineticRxns_.begin();
       i != generalKineticRxns_.end(); i++) {
    i->update_rates(primarySpecies_);
  }

  // add mineral saturation and rate calculations here
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->Update(primarySpecies_, minerals_);
  }
  // add multirate kinetic surface complexation reaction quotient calculations
  // here

}  // end updateKineticChemistry()

void Beaker::addKineticChemistryToResidual(std::vector<double> *residual) {
  // loop over general kinetic reactions and add rates
  for (std::vector<GeneralRxn>::iterator i = generalKineticRxns_.begin();
       i != generalKineticRxns_.end(); i++) {
    i->addContributionToResidual(residual, por_sat_den_vol());
  }

  // add mineral mineral contribution to residual here.  units = mol/sec.
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->AddContributionToResidual(minerals_, por_sat_den_vol(), residual);
  }

  // add multirate kinetic surface complexation contribution to residual here.

}  // end addKineticChemistryToResidual()

void Beaker::addKineticChemistryToJacobian(Block* J) {
  // loop over general kinetic reactions and add rates
  for (std::vector<GeneralRxn>::iterator i = generalKineticRxns_.begin();
       i != generalKineticRxns_.end(); i++) {
    i->addContributionToJacobian(J, primarySpecies_, por_sat_den_vol());
  }

  // add mineral mineral contribution to Jacobian here.  units = kg water/sec.
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->AddContributionToJacobian(primarySpecies_, minerals_, por_sat_den_vol(), J);
  }

  // add multirate kinetic surface complexation contribution to Jacobian here.

}  // end addKineticChemistryToJacobian()



void Beaker::addAccumulation(std::vector<double> *residual) {
  addAccumulation(total_, total_sorbed_, residual);
}

void Beaker::addAccumulation(const std::vector<double>& total,
                             const std::vector<double>& total_sorbed,
                             std::vector<double> *residual) {
  // accumulation_coef = porosity*saturation*volume*1000./dt
  // units = (mol solute/L water)*(m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  //         (m^3 bulk)*(1000L water/m^3 water)/(sec) = mol/sec
  // 1000.d0 converts vol from m^3 -> L
  // all residual entries should be in mol/sec
  for (unsigned int i = 0; i < total.size(); i++) {
    (*residual)[i] += aqueous_accumulation_coef() * total[i];
  }

  // add accumulation term for equilibrium sorption (e.g. Kd, surface
  // complexation) here
  // accumulation_coef = porosity*saturation*volume*1000./dt
  // units = (mol solute/L water)*(m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  //         (m^3 bulk)*(1000L water/m^3 water)/(sec) = mol/sec
  // 1000.d0 converts vol from m^3 -> L
  // all residual entries should be in mol/sec
  for (unsigned int i = 0; i < total_sorbed.size(); i++) {
    (*residual)[i] += sorbed_accumulation_coef() * total_sorbed[i];
  }

}  // end calculateAccumulation()

void Beaker::addAccumulationDerivative(Block* J) {
  addAccumulationDerivative(J, dtotal_, dtotal_sorbed_);
}  // end calculateAccumulationDerivative()

void Beaker::addAccumulationDerivative(Block* J,
                                       Block* dtotal,
                                       Block* dtotal_sorbed) {
  // accumulation_coef = porosity*saturation*volume*1000./dt
  // units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*(m^3 bulk)/(sec)
  //         *(kg water/L water)*(1000L water/m^3 water) = kg water/sec
  // all Jacobian entries should be in kg water/sec
  J->addValues(dtotal, aqueous_accumulation_coef());

  // add accumulation derivative term for equilibrium sorption
  // (e.g. Kd, surface complexation) here
  J->addValues(dtotal_sorbed, sorbed_accumulation_coef());

}  // end calculateAccumulationDerivative()

void Beaker::calculateFixedAccumulation(const std::vector<double>& total,
                                        const std::vector<double>& total_sorbed,
                                        std::vector<double> *fixed_accumulation) {
  for (unsigned int i = 0; i < total.size(); i++) {
    (*fixed_accumulation)[i] = 0.;
  }
  addAccumulation(total, total_sorbed, fixed_accumulation);
}  // end calculateAccumulation()

void Beaker::calculateResidual(std::vector<double> *residual,
                               const std::vector<double>& fixed_accumulation) {
  status_.num_rhs_evaluations++;
  // subtract fixed porition
  for (int i = 0; i < ncomp(); i++) {
    (*residual)[i] = -fixed_accumulation[i];
  }

  // accumulation adds in equilibrium chemistry
  addAccumulation(residual);

  // kinetic reaction contribution to residual
  addKineticChemistryToResidual(residual);


}  // end calculateResidual()

void Beaker::calculateJacobian(Block* J) {
  status_.num_jacobian_evaluations++;
  // must calculate derivatives with
  calculateDTotal();

  // zero Jacobian
  J->zero();
  // add in derivatives for equilibrium chemistry
  addAccumulationDerivative(J);

  // add in derivatives for kinetic chemistry
  addKineticChemistryToJacobian(J);

}  // end calculateJacobian()

void Beaker::scaleRHSAndJacobian(double* rhs, Block* J) {

  for (int i = 0; i < J->getSize(); i++) {
    double max = J->getRowAbsMax(i);
    if (max > 1.) {
      double scale = 1. / max;
      rhs[i] *= scale;
      J->scaleRow(i, scale);
    }
  }
}  // end scaleRHSAndJacobian()

void Beaker::scaleRHSAndJacobian(std::vector<double> &rhs, Block* J) {

  for (int i = 0; i < J->getSize(); i++) {
    double max = J->getRowAbsMax(i);
    if (max > 1.) {
      double scale = 1. / max;
      rhs[i] *= scale;
      J->scaleRow(i, scale);
    }
  }
}  // end scaleRHSAndJacobian()

void Beaker::updateMolalitiesWithTruncation(std::vector<double> &update,
                                            std::vector<double> &prev_solution,
                                            double max_change) {
  // truncate the update to max_change
  for (int i = 0; i < ncomp(); i++) {
    if (update[i] > max_change) {
      update[i] = max_change;
    } else if (update[i] < -max_change) {
      update[i] = -max_change;
    }
    // store the previous solution
    prev_solution[i] = primarySpecies_[i].molality();
    // update primary species molalities (log formulation)
    double molality = prev_solution[i] * std::exp(-update[i]);
    primarySpecies_[i].update(molality);
  }
}  // end updateMolalitiesWithTruncation()

double Beaker::calculateMaxRelChangeInMolality(const std::vector<double>& prev_molal) {
  double max_rel_change = 0.0;
  for (int i = 0; i < ncomp(); i++) {
    double delta = fabs(primarySpecies_[i].molality() - prev_molal[i]) / prev_molal[i];
    max_rel_change = delta > max_rel_change ? delta : max_rel_change;
  }
  return max_rel_change;
}  // end calculateMaxRelChangeInMolality()

void Beaker::solveLinearSystem(Block* A, std::vector<double> &b) {
  // LU direct solve
#ifndef GLENN
  double D;
  // allocate pivoting array for LU
  int* indices = new int[ncomp()];
  /* TODO:this is a memory leak? we are new'ing a bunch of memory
     but not freeing it.  need to change the object level indices
     variable into an int*, then allocate correctly in resize() and
     free in the destructor */
  ludcmp(A->getValues(), ncomp(), indices, &D);
  lubksb(A->getValues(), ncomp(), indices, b);
  delete[] indices;
#else
  solver->Solve(A, b);
#endif
}  // end solveLinearSystem

void Beaker::CheckChargeBalance(const std::vector<double>& aqueous_totals) {
  double charge_balance = 0.0;
  for (unsigned int i = 0; i < aqueous_totals.size(); i++) {
    charge_balance += aqueous_totals.at(i) * primarySpecies_.at(i).charge();
  }
  if (charge_balance > tolerance()) {
    if (verbosity() > kTerse) {
      std::cout << "WARNING: Beaker::CheckChargeBalance() :\n"
                << "         charge balance = " << std::scientific
                << charge_balance << std::fixed << std::endl;
    }
  }
}  // end CheckChargeBalance()

int Beaker::ReactionStep(Beaker::BeakerComponents* components,
                         const Beaker::BeakerParameters& parameters,
                         double dt) {
  /*
  ** Note: the parameter components is modified by this function.
  ** initially it contains the initial component concentrations.
  ** on return it contains the modified values of the components.
  */

  // update class paramters
  // water_density [kg/m^3]
  // volume [m^3]
  ResetStatus();
  updateParameters(parameters, dt);
  CheckChargeBalance(components->total);

  // initialize free-ion concentrations, these are actual
  // concentrations, so the value must be positive or ln(free_ion)
  // will return a nan!
  if (components->free_ion.size() > 0) {
    initializeMolalities(components->free_ion);
    // testing only    initializeMolalities(1.e-9);
  }

  // store current molalities
  for (int i = 0; i < ncomp(); i++) {
    prev_molal[i] = primarySpecies_[i].molality();
  }

  for (unsigned int m = 0; m < minerals_.size(); m++) {
    minerals_.at(m).set_volume_fraction(components->minerals.at(m));
    minerals_.at(m).UpdateSurfaceAreaFromVolumeFraction(volume());
  }

  // initialize to a large number (not necessary, but safe)
  double max_rel_change = 1.e20;
  // iteration counter
  unsigned int num_iterations = 0;

  // lagging activity coefficients by a time step in this case
  updateActivityCoefficients();

  // calculate portion of residual at time level t
  calculateFixedAccumulation(components->total, components->total_sorbed,
                             &fixed_accumulation);

  do {

    // update equilibrium and kinetic chemistry (rates, ion activity
    // products, etc.)
    updateEquilibriumChemistry();
    updateKineticChemistry();

    // units of residual: mol/sec
    calculateResidual(&residual, fixed_accumulation);
    // units of Jacobian: kg water/sec
    calculateJacobian(J);
    // therefore, units of solution: mol/kg water (change in molality)

    for (int i = 0; i < ncomp(); i++) {
      rhs[i] = residual[i];
    }

    if (verbosity() == kDebugBeaker) {
      print_linear_system("before scale", J, rhs);
    }

    // scale the Jacobian
    scaleRHSAndJacobian(rhs, J);

    if (verbosity() == kDebugBeaker) {
      print_linear_system("after scale", J, rhs);
    }
    // for derivatives with respect to ln concentration, scale columns
    // by primary species concentrations
    for (int i = 0; i < ncomp(); i++) {
      J->scaleColumn(i, primarySpecies_[i].molality());
    }

    if (verbosity() == kDebugBeaker) {
      print_linear_system("before solve", J, rhs);
    }

    // solve J dlnc = r
    solveLinearSystem(J, rhs);

    if (verbosity() == kDebugBeaker) {
      print_linear_system("after solve", NULL, rhs);
    }

    // units of solution: mol/kg water (change in molality)
    // calculate update truncating at a maximum of 5 in nat log space
    // update with exp(-dlnc)
    updateMolalitiesWithTruncation(rhs, prev_molal, 5.);
    // calculate maximum relative change in concentration over all species
    max_rel_change = calculateMaxRelChangeInMolality(prev_molal);

    if (verbosity() == kDebugBeaker) {
      for (int i = 0; i < ncomp(); i++)
        std::cout << primarySpecies_[i].name() << " " <<
            primarySpecies_[i].molality() << " " << total_[i] << "\n";
    }

    num_iterations++;

    // geh
    //    if (num_iterations >= 100) {
    //      for (int i = 0; i < ncomp(); i++)
    //        std::cout << primarySpecies_[i].name() << " " <<
    //                  primarySpecies_[i].molality() << " " << total_[i] << "\n";
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
  updateEquilibriumChemistry();
  UpdateComponents(components);
  ValidateSolution();

  return num_iterations;

}  // end ReactionStep()

void Beaker::ValidateSolution() {
  // TODO: what checks can we to to verify that the current solution is good?
  // TODO: check for nan or inf's as a sign that the step was too big?

  // TODO: negative total's (H+) are OK...

  // charge balance is error or warning...?
  CheckChargeBalance(total_);

  // negative mineral volume fractions are bad...
  for (unsigned int m = 0; m < minerals_.size(); m++) {
    if (minerals_.at(m).volume_fraction() < 0) {
      std::ostringstream error_stream;
      error_stream << "Beaker::ValidateSolution(): \n";
      error_stream << "   mineral " << minerals_.at(m).name()
                   << " volume_fraction is negative." << std::endl;
      Exceptions::amanzi_throw(ChemistryInvalidSolution(error_stream.str()));
    }
  }

}  // end ValidateSolution()


// if no water density provided, default is 1000.0 kg/m^3
int Beaker::Speciate(const Beaker::BeakerComponents& components,
                     const Beaker::BeakerParameters& parameters) {
  double speciation_tolerance = 1.e-12;
  double residual_tolerance = 1.e-12;
  ResetStatus();
  updateParameters(parameters, 0.0);
  CheckChargeBalance(components.total);
  // initialize free-ion concentrations, these are actual
  // concentrations, so the value must be positive or ln(free_ion)
  // will return a nan!
  if (components.free_ion.size() > 0) {
    initializeMolalities(components.free_ion);
  } else {
    initializeMolalities(1.e-9);
  }

  for (unsigned int m = 0; m < minerals_.size(); m++) {
    minerals_.at(m).set_volume_fraction(components.minerals.at(m));
    minerals_.at(m).UpdateSurfaceAreaFromVolumeFraction(volume());
  }

  // store current molalities
  for (int i = 0; i < ncomp(); i++) {
    prev_molal[i] = primarySpecies_[i].molality();
  }

  double max_rel_change;
  double max_residual;
  unsigned int num_iterations = 0;
  bool calculate_activity_coefs = false;

  do {

    updateActivityCoefficients();
    updateEquilibriumChemistry();
    calculateDTotal();

    // calculate residual
    // units of residual: mol/sec
    for (int i = 0; i < ncomp(); i++) {
      residual[i] = total_[i] - components.total[i];
    }
    // add derivatives of total with respect to free to Jacobian
    // units of Jacobian: kg water/sec
    J->zero();
    calculateDTotal();
    J->addValues(0, 0, dtotal_);

    for (int i = 0; i < ncomp(); i++) {
      rhs[i] = residual[i];
    }

    // geh
    DisplayResults();

    if (verbosity() == kDebugBeaker) {
      print_linear_system("before scale", J, rhs);
    }

    // scale the Jacobian
    scaleRHSAndJacobian(rhs, J);

    if (verbosity() == kDebugBeaker) {
      print_linear_system("after scale", J, rhs);
    }
    // for derivatives with respect to ln concentration, scale columns
    // by primary species concentrations
    for (int i = 0; i < ncomp(); i++) {
      J->scaleColumn(i, primarySpecies_[i].molality());
    }

    if (verbosity() == kDebugBeaker) {
      print_linear_system("before solve", J, rhs);
    }

    // call solver
    solveLinearSystem(J, rhs);

    if (verbosity() == kDebugBeaker) {
      print_linear_system("after solve", NULL, rhs);
    }

    // calculate update truncating at a maximum of 5 in log space
    updateMolalitiesWithTruncation(rhs, prev_molal, 5.);
    // calculate maximum relative change in concentration over all species
    max_rel_change = calculateMaxRelChangeInMolality(prev_molal);

    if (verbosity() == kDebugBeaker) {
      for (int i = 0; i < ncomp(); i++) {
        std::cout << primarySpecies_[i].name() << " "
                  << primarySpecies_[i].molality() << " " << total_[i] << "\n";
      }
    }

    num_iterations++;

    max_residual = 0.;
    for (int i = 0; i < ncomp(); i++)
      if (fabs(residual[i]) > max_residual) {
        max_residual = fabs(residual[i]);
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
  updateEquilibriumChemistry();
  status_.num_newton_iterations = num_iterations;
  if (max_rel_change < tolerance()) {
    status_.converged = true;
  }

  if (verbosity() >= kDebugBeaker) {
    std::cout << "Beaker::speciate: status.num_rhs_evaluations: " << status_.num_rhs_evaluations << std::endl;
    std::cout << "Beaker::speciate: status.num_jacobian_evaluations: " << status_.num_jacobian_evaluations << std::endl;
    std::cout << "Beaker::speciate: status.num_newton_iterations: " << status_.num_newton_iterations << std::endl;
    std::cout << "Beaker::speciate: status.converged: " << status_.converged << std::endl;
  }
  return num_iterations;
}  // end Speciate()


void Beaker::UpdateComponents(Beaker::BeakerComponents* components) {
  for (int i = 0; i < ncomp(); i++) {
    components->total[i] = total_[i];
    if (components->free_ion.size() > 0) {
      components->free_ion[i] = primarySpecies_[i].molality();
    }
    if (components->total_sorbed.size() > 0 &&
        total_sorbed_.size() > 0) {
      components->total_sorbed[i] = total_sorbed_[i];
    }
  }
  for (int m = 0; m < minerals_.size(); m++) {
    components->minerals[m] = minerals_.at(m).volume_fraction();
  }
  for (int i = 0; i < ion_exchange_sites_.size(); i++) {
    components->ion_exchange_sites[i] = ion_exchange_sites_.at(i).cation_exchange_capacity();
  }
}  // end UpdateComponents()

void Beaker::CopyComponents(const Beaker::BeakerComponents& from,
                            Beaker::BeakerComponents* to) {
  if (to->total.size() != from.total.size()) {
    to->total.resize(from.total.size());
  }
  if (to->free_ion.size() != from.free_ion.size()) {
    to->free_ion.resize(from.free_ion.size());
  }
  if (to->total_sorbed.size() != from.total_sorbed.size()) {
    to->total_sorbed.resize(from.total_sorbed.size());
  }
  if (to->minerals.size() != from.minerals.size()) {
    to->minerals.resize(from.minerals.size());
  }
  if (to->ion_exchange_sites.size() != from.ion_exchange_sites.size()) {
    to->ion_exchange_sites.resize(from.ion_exchange_sites.size());
  }
  for (unsigned int i = 0; i < from.total.size(); i++) {
    to->total[i] = from.total[i];
  }
  for (unsigned int i = 0; i < from.free_ion.size(); i++) {
    to->free_ion[i] = from.free_ion[i];
  }
  for (unsigned int i = 0; i < from.total_sorbed.size(); i++) {
    to->total_sorbed[i] = from.total_sorbed[i];
  }
  for (unsigned int i = 0; i < from.minerals.size(); i++) {
    to->minerals[i] = from.minerals[i];
  }
  for (unsigned int i = 0; i < from.ion_exchange_sites.size(); i++) {
    to->ion_exchange_sites[i] = from.ion_exchange_sites[i];
  }
}  // end Beaker::CopyComponents()

/*******************************************************************************
 **
 **  Output related functions
 **
 *******************************************************************************/
void Beaker::GetPrimaryNames(std::vector<std::string>* names) const {
  names->clear();
  for (std::vector<Species>::const_iterator primary = primarySpecies_.begin();
       primary != primarySpecies_.end(); primary++) {
    names->push_back(primary->name());
  }
}  // end GetPrimaryNames()

int Beaker::GetPrimaryIndex(const std::string& name) const {
  int index = -1;
  for (std::vector<Species>::const_iterator primary = primarySpecies_.begin();
       primary != primarySpecies_.end(); primary++) {
    if (primary->name() == name) {
      index = primary->identifier();
    }
  }
  return index;
}  // end GetPrimaryIndex()

void Beaker::display(void) const {
  std::cout << "----- Beaker description ------" << std::endl;
  std::cout << "Primary Species:" << std::endl;
  for (std::vector<Species>::const_iterator primary = primarySpecies_.begin();
       primary != primarySpecies_.end(); primary++) {
    primary->display();
  }
  std::cout << std::endl;
  std::cout << "Aqueous Equilibrium Complexes:" << std::endl;
  for (std::vector<AqueousEquilibriumComplex>::const_iterator aec = aqComplexRxns_.begin();
       aec != aqComplexRxns_.end(); aec++) {
    aec->display();
  }
  std::cout << "-------------------------------------" << std::endl;
}  // end display()

void Beaker::Display(void) const {
  std::cout << "-- Beaker description ------------------------------------------------"
            << std::endl;
  if (verbosity() >= kVerbose) {
    DisplayParameters();
  }

  DisplayPrimary();

  DisplayAqueousEquilibriumComplexes();

  DisplayMinerals();

  DisplayMineralKinetics();

  DisplayIonExchangeSites();

  DisplayIonExchangeComplexes();

  DisplaySurfaceSites();

  DisplaySurfaceComplexes();

  std::cout << "----------------------------------------------------------------------"
            << std::endl;

}  // end Display()

void Beaker::DisplayParameters(void) const {
  // units....
  std::cout << "---- Parameters" << std::endl;
  // std::cout << "    thermo_database_file: " << thermo_database_file << std::endl;
  std::cout << "    tolerance: " << tolerance() << std::endl;
  std::cout << "    max_iterations :" << max_iterations() << std::endl;

  std::cout << "    activity model: " << activity_model_->name() << std::endl;

  std::cout << "    porosity: " << porosity() << " [-]" << std::endl;
  std::cout << "    water saturation: " << saturation() << " [-]" << std::endl;
  std::cout << "    water density: " << water_density_kg_m3() << " [kg m^-3]" << std::endl;
  std::cout << "    volume: " << volume() << " [m^3]" << std::endl;
  std::cout << std::endl;

}  // end DisplayParameters()

void Beaker::DisplayPrimary(void) const {
  std::cout << "---- Primary Species" << std::endl;
  std::cout << std::setw(15) << "Species"
            << std::setw(10) << "Charge"
            << std::setw(10) << "GMW"
            << std::setw(10) << "D-H a0"
            << std::endl;
  for (std::vector<Species>::const_iterator primary = primarySpecies_.begin();
       primary != primarySpecies_.end(); primary++) {
    primary->Display();
  }
  std::cout << std::endl;
}  // end DisplayPrimary()

void Beaker::DisplayAqueousEquilibriumComplexes(void) const {
  std::cout << "---- Aqueous Equilibrium Complexes" << std::endl;
  std::cout << std::setw(12) << "Reaction"
            << std::setw(38) << "log Keq"
            << std::setw(8) << "Charge"
            << std::setw(10) << "GMW"
            << std::setw(8) << "D-H a0"
            << std::endl;
  for (std::vector<AqueousEquilibriumComplex>::const_iterator aec = aqComplexRxns_.begin();
       aec != aqComplexRxns_.end(); aec++) {
    aec->Display();
  }
  std::cout << std::endl;
}  // end DisplayAqueousEquilibriumComplexes()

void Beaker::DisplayMinerals(void) const {
  if (minerals_.size() > 0) {
    std::cout << "---- Minerals" << std::endl;
    std::cout << std::setw(12) << "Reaction"
              << std::setw(38) << "log_Keq"
              << std::setw(13) << "molar volume"
              << std::setw(13) << "GMW"
              << std::setw(13) << "SSA"
              << std::setw(13) << "Area"
              << std::setw(13) << "Vfrac"
              << std::endl;
    std::cout << std::setw(12) << " "
              << std::setw(38) << " "
              << std::setw(13) << "[cm^3/mol]"
              << std::setw(13) << "[g/mol]"
              << std::setw(13) << "[m^2/g]"
              << std::setw(13) << "[m^2]"
              << std::setw(13) << "[-]"
              << std::endl;
    for (std::vector<Mineral>::const_iterator m = minerals_.begin();
         m != minerals_.end(); m++) {
      m->Display();
    }
    std::cout << std::endl;
  }
}  // end DisplayMinerals()

void Beaker::DisplayMineralKinetics(void) const {
  if (mineral_rates_.size() > 0) {
    std::cout << "---- Mineral Kinetics" << std::endl;
    for (std::vector<KineticRate*>::const_iterator m = mineral_rates_.begin();
         m != mineral_rates_.end(); m++) {
      (*m)->Display();
    }
    std::cout << std::endl;
  }
}  // end DisplayMineralKinetics()

void Beaker::DisplayIonExchangeSites(void) const {
  if (ion_exchange_sites_.size() > 0) {
    std::cout << "---- Ion Exchange Sites" << std::endl;
    std::cout << std::setw(15) << "Species"
              << std::setw(15) << "Location"
              << std::setw(10) << "Charge"
              << std::setw(10) << "CEC"
              << std::endl;
    std::vector<IonExchangeSite>::const_iterator exchanger;
    for (exchanger = ion_exchange_sites_.begin();
         exchanger != ion_exchange_sites_.end(); exchanger++) {
      exchanger->Display();
    }
    std::cout << std::endl;
  }
}  // end DisplayIonExchangeSites()

void Beaker::DisplayIonExchangeComplexes(void) const {
  if (ion_exchange_rxns_.size() > 0) {
    std::cout << "---- Ion Exchange Complexes" << std::endl;
    std::cout << std::setw(12) << "Reaction"
              << std::setw(38) << "log Keq"
              << std::endl;
    std::vector<IonExchangeComplex>::const_iterator iec;
    for (iec = ion_exchange_rxns_.begin();
         iec != ion_exchange_rxns_.end(); iec++) {
      iec->Display();
    }
    std::cout << std::endl;
  }
}  // end DisplayIonExchangeComplexes()

void Beaker::DisplaySurfaceSites(void) const {
  if (total_sorbed_.size() > 0) {
    std::cout << "---- Surface Sites" << std::endl;
    std::cout << std::setw(15) << "Species"
              << std::setw(15) << "Site Density"
              << std::endl;
    std::vector<SurfaceComplexationRxn>::const_iterator s;
    for (s = surfaceComplexationRxns_.begin();
         s != surfaceComplexationRxns_.end(); s++) {
      s->DisplaySite();
    }
    std::cout << std::endl;
  }
}  // end DisplaySurfaceSites()

void Beaker::DisplaySurfaceComplexes(void) const {
  if (total_sorbed_.size() > 0) {
    std::cout << "---- Surface Complexes" << std::endl;
    std::cout << std::setw(12) << "Reaction"
              << std::setw(38) << "log Keq"
              << std::setw(10) << "charge"
              << std::endl;
    std::vector<SurfaceComplexationRxn>::const_iterator s;
    for (s = surfaceComplexationRxns_.begin();
         s != surfaceComplexationRxns_.end(); s++) {
      s->DisplayComplexes();
    }
    std::cout << std::endl;
  }
}  // end DisplaySurfaceComplexes()

void Beaker::DisplayComponents(const Beaker::BeakerComponents& components) const {
  std::cout << "--- Input Components -------------------------------------------------"
            << std::endl;
  std::cout << "---- Aqueous Components" << std::endl;
  std::cout << std::setw(15) << "Name"
            << std::setw(15) << "Molarity"
            << std::setw(15) << "Molality"
      // << std::setw(15) << "Free Ion" // TODO: uncomment and update test results
            << std::endl;
  for (int i = 0; i < ncomp(); i++) {
    std::cout << std::setw(15) << primarySpecies_.at(i).name()
              << std::scientific << std::setprecision(5)
              << std::setw(15) << components.total.at(i) / water_density_kg_L()
              << std::setw(15) << components.total.at(i)
        // << std::setw(15) << components.free_ion.at(i) // TODO: uncomment and update test results
              << std::endl;
  }

  if (minerals_.size() > 0) {
    std::cout << "---- Mineral Components" << std::endl;
    std::cout << std::setw(15) << "Name"
              << std::setw(15) << "Vol. frac" << std::endl;
    for (unsigned int m = 0; m < minerals_.size(); m++) {
      std::cout << std::setw(15) << minerals_.at(m).name()
                << std::setw(15) << std::fixed << std::setprecision(5)
                << components.minerals.at(m) << std::endl;
    }
  }

  if (total_sorbed_.size() > 0) {
    std::cout << "---- Surface Complex Sites" << std::endl;
    std::cout << std::setw(15) << "Name"
              << std::setw(15) << "Molarity" << std::endl;
    for (int i = 0; i < ncomp(); i++) {
      std::cout << std::setw(15) << primarySpecies_.at(i).name()
                << std::scientific << std::setprecision(5)
                << std::setw(15) << components.total_sorbed.at(i)
                << std::endl;
    }
  }
  std::cout << "----------------------------------------------------------------------"
            << std::endl;

}  // end DisplayComponents

void Beaker::DisplayResults(void) const {
  std::cout << std::endl;
  std::cout << "-- Solution ----------------------------------------------------------"
            << std::endl;
  std::cout << "---- Components " << std::endl;
  std::cout << std::setw(15) << "Name"
            << std::setw(15) << "Molarity"
            << std::setw(15) << "Molality"
            << std::endl;
  for (int i = 0; i < ncomp(); i++) {
    std::cout << std::setw(15) << primarySpecies_.at(i).name()
              << std::scientific << std::setprecision(5)
              << std::setw(15) << total_[i] / water_density_kg_L()
              << std::setw(15) << total_[i]
              << std::endl;
  }

  std::cout << "---- Change Balance " << std::endl;
  double charge_balance_molal = 0.0;
  for (int i = 0; i < ncomp(); i++) {
    charge_balance_molal += primarySpecies_.at(i).charge() * total_[i];
  }
  std::cout << std::setw(15) << " "
            << std::scientific << std::setprecision(5)
            << std::setw(15) << " "
            << std::setw(15) << charge_balance_molal
            << std::endl;

  std::cout << "---- Species " << std::endl;

  primarySpecies_.at(0).DisplayResultsHeader();
  for (int i = 0; i < ncomp(); i++) {
    primarySpecies_.at(i).DisplayResults();
  }

  // same header info as primaries....
  for (int i = 0; i < (int)aqComplexRxns_.size(); i++) {
    aqComplexRxns_.at(i).DisplayResults();
  }

  if (minerals_.size() > 0) {
    std::cout << "---- Minerals " << std::endl;
    minerals_[0].DisplayResultsHeader();
    for (unsigned int i = 0; i < minerals_.size(); i++) {
      minerals_[i].DisplayResults();
    }
  }

  if (ion_exchange_sites_.size() > 0) {
    std::cout << "---- Ion Exchange Sites " << std::endl;
    ion_exchange_sites_[0].DisplayResultsHeader();
    for (unsigned int i = 0; i < ion_exchange_sites_.size(); i++) {
      ion_exchange_sites_[i].DisplayResults();
    }
  }

  if (ion_exchange_rxns_.size() > 0) {
    std::cout << "---- Ion Exchange Complexes " << std::endl;
    ion_exchange_rxns_[0].DisplayResultsHeader();
    for (unsigned int i = 0; i < ion_exchange_rxns_.size(); i++) {
      ion_exchange_rxns_[i].DisplayResults();
    }
  }

  if (surfaceComplexationRxns_.size() > 0) {
    std::cout << "---- Surface Complexation Reactions " << std::endl;
    for (unsigned int i = 0; i < surfaceComplexationRxns_.size(); i++) {
      surfaceComplexationRxns_[i].DisplayResultsHeader();
      surfaceComplexationRxns_[i].DisplayResults();
    }
  }

  std::cout << "----------------------------------------------------------------------"
            << std::endl << std::endl;
}  // end DisplayResults()

void Beaker::DisplayTotalColumnHeaders(void) const {
  std::cout << std::setw(15) << "Time (s)";
  for (int i = 0; i < ncomp(); i++) {
    std::cout << std::setw(15) << primarySpecies_.at(i).name();
  }
  std::cout << std::endl;
}  // end DisplayTotalColumnHeaders()

void Beaker::DisplayTotalColumns(const double time, const std::vector<double>& total) const {
  std::cout << std::scientific << std::setprecision(5) << std::setw(15);
  std::cout << time;
  for (int i = 0; i < ncomp(); i++) {
    std::cout << std::setw(15) << total.at(i);
  }
  std::cout << std::endl;

}  // end DisplayTotalColumns()


void Beaker::print_results(void) const {
  // output for testing purposes
  std::cout << std::endl;
  std::cout << "----- Solution ----------------------" << std::endl;
  std::cout << "Primary Species ---------------------\n";
  for (int i = 0; i < ncomp(); i++) {
    std::cout << "   " << primarySpecies_.at(i).name() << std::endl;
    std::cout << "        Total: " << total_[i] << std::endl;
    std::cout << "     Free-Ion: " << primarySpecies_.at(i).molality() << std::endl;
    std::cout << "Activity Coef: " << primarySpecies_.at(i).act_coef() << std::endl;
    std::cout << "     Activity: " << primarySpecies_.at(i).activity() << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Secondary Species -------------------\n";
  for (int i = 0; i < (int)aqComplexRxns_.size(); i++) {
    std::cout << "   " << aqComplexRxns_.at(i).name() << std::endl;
    std::cout << "     Free-Ion: " << aqComplexRxns_.at(i).molality() << std::endl;
    std::cout << "Activity Coef: " << aqComplexRxns_.at(i).act_coef() << std::endl;
    std::cout << "     Activity: " << aqComplexRxns_.at(i).activity() << std::endl;
  }
  std::cout << "-------------------------------------\n";
  std::cout << std::endl;

}  // end print_results()

void Beaker::print_results(double time) const {
  if (time < 1.e-40) {
    std::cout << "Time\t";
    for (int i = 0; i < ncomp(); i++) {
      std::cout << primarySpecies_.at(i).name() << " (total)\t";
      std::cout << primarySpecies_.at(i).name() << " (free-ion)\t";
    }
    std::cout << std::endl;
  }
  // output for testing purposes
  std::cout << time << "\t";
  for (int i = 0; i < ncomp(); i++) {
    std::cout << total_[i] << "\t";
    std::cout << primarySpecies_.at(i).molality() << "\t";
  }
  std::cout << std::endl;
}  // end print_results()

void Beaker::print_linear_system(const std::string& s, Block* A,
                                 std::vector<double> vector) {
  std::cout << s << std::endl;
  for (int i = 0; i < (int)vector.size(); i++) {
    std::cout << "RHS: " << primarySpecies_.at(i).name() << " " << vector[i] << std::endl;
  }
  if (A) {
    A->print();
  }
}  // end print_linear_system()
