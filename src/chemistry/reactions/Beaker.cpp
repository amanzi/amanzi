#include <cstdlib>

#include "Beaker.hpp"
#include "ActivityModelFactory.hpp"
#include "KineticRate.hpp"
#include "MineralKineticsFactory.hpp"
#include "Verbosity.hpp"

// solver defaults
const double Beaker::tolerance_default = 1.0e-12;
const unsigned int Beaker::max_iterations_default = 250;
// default physical parameters
const double Beaker::porosity_default = 1.0; // [-]
const double Beaker::saturation_default = 1.0; // [-]
const double Beaker::water_density_kg_m3_default = 1000.0; 
const double Beaker::volume_default = 1.0; // [m^3]

Beaker::Beaker() 
  : verbosity_(kSilent),
    tolerance_(tolerance_default),
    max_iterations_(max_iterations_default),
    ncomp_(0),
    dtotal_(NULL),
    porosity_(porosity_default),
    saturation_(saturation_default),
    water_density_kg_m3_(water_density_kg_m3_default),
    water_density_kg_L_(1.0),
    volume_(volume_default),
    dt_(0.0),
    accumulation_coef_(0.0), 
    por_sat_den_vol_(0.0),
    activity_model_(NULL),
    J(NULL)
{
  aqComplexRxns_.clear();
  generalKineticRxns_.clear();
  mineral_rates_.clear();
} // end Beaker() constructor

Beaker::~Beaker() 
{
  if (dtotal_ != NULL) {
    delete dtotal_;
  }

  if (J != NULL) {
    delete J;
  }

  if (activity_model_ != NULL) {
    delete activity_model_;
  }

  if (mineral_rates_.size() != 0) {
    for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
         rate != mineral_rates_.end(); rate++) {
      if ((*rate) != NULL) {
        delete (*rate);
      }
    }    
  }

} // end Beaker destructor

void Beaker::resize() {
  if (dtotal_) delete dtotal_;
  dtotal_ = new Block(ncomp());
  fixed_accumulation.resize(ncomp());
  residual.resize(ncomp());
  prev_molal.resize(ncomp());
  total_.resize(ncomp());

  if (J) delete J;
  J = new Block(ncomp());
  rhs.resize(ncomp());
  indices.resize(ncomp());

} // end resize()

void Beaker::resize(int n) {
  ncomp(n);
  resize();
} // end resize()

void Beaker::setup(std::vector<double> &total, 
                   const std::string mineral_kinetics_file)
{
  static_cast<void>(total);
  static_cast<void>(mineral_kinetics_file);
} // end setup()

void Beaker::SetupActivityModel(std::string model)
{
  if (model != ActivityModelFactory::unit &&
      model != ActivityModelFactory::debye_huckel) {
    model = ActivityModelFactory::unit;
  }
  if (activity_model_ != NULL) {
    delete activity_model_;
  }
  ActivityModelFactory amf;
  
  activity_model_ = amf.Create(model);
  if (verbosity() >= kVerbose) {
    activity_model_->Display();
  }
}  // end SetupActivityModel() 

void Beaker::SetupMineralKinetics(const std::string mineral_kinetics_file)
{
  if (mineral_kinetics_file.size()) {
    MineralKineticsFactory mineral_kinetics_factory;
    mineral_kinetics_factory.verbosity(verbosity());
    mineral_rates_ = mineral_kinetics_factory.Create(mineral_kinetics_file, 
                                                     primarySpecies_);
    
    if (verbosity() >= kVerbose) {
      for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
           rate != mineral_rates_.end(); rate++) {
        (*rate)->Display();
      }
    }
  }

}  // end SetupMineralKinetics()

void Beaker::addPrimarySpecies(Species s) 
{
  primarySpecies_.push_back(s);
} // end addPrimarySpecies()

void Beaker::addAqueousEquilibriumComplex(AqueousEquilibriumComplex c) 
{
  aqComplexRxns_.push_back(c);
} // end addAqueousEquilibriumComplex()

void Beaker::addGeneralRxn(GeneralRxn r) 
{
  generalKineticRxns_.push_back(r);
} // end addGeneralRxn()

Beaker::BeakerParameters Beaker::GetDefaultParameters(void)
{
  Beaker::BeakerParameters parameters;

  parameters.tolerance = tolerance_default;
  parameters.max_iterations = max_iterations_default;
 
  parameters.porosity = porosity_default;
  parameters.saturation = saturation_default;
  parameters.water_density = water_density_kg_m3_default; // kg / m^3
  parameters.volume = volume_default; // m^3
 
  return parameters;
} // end GetDefaultParameters()
    
Beaker::BeakerParameters Beaker::GetCurrentParameters(void)
{
  Beaker::BeakerParameters parameters;

  parameters.tolerance = tolerance();
  parameters.max_iterations = max_iterations();
 
  parameters.porosity = porosity();
  parameters.saturation = saturation();
  parameters.water_density = water_density_kg_m3(); // kg / m^3
  parameters.volume = volume(); // m^3
 
  return parameters;
} // end GetCurrentParameters()
    

void Beaker::updateParameters(const Beaker::BeakerParameters& parameters, 
                              double delta_t)
{
  tolerance(parameters.tolerance);
  max_iterations(parameters.max_iterations);
  porosity(parameters.porosity);
  water_density_kg_m3(parameters.water_density); // den = [kg/m^3]
  saturation(parameters.saturation);
  volume(parameters.volume); // vol = [m^3]
  dt(delta_t); // delta time = [sec]
  update_accumulation_coef();
  update_por_sat_den_vol();
} // end updateParameters()

void Beaker::update_accumulation_coef(void)
{
  accumulation_coef(porosity() * saturation() * volume() * 1000.0 / dt());
} // end update_accumulation_coef

void Beaker::update_por_sat_den_vol(void)
{
  por_sat_den_vol(porosity() * saturation() * water_density_kg_m3() * volume()); 
} // end update_por_sat_den_vol()


void Beaker::updateActivityCoefficients() {

  //return;
  activity_model_->CalculateIonicStrength(primarySpecies_,
                                          aqComplexRxns_);
  activity_model_->CalculateActivityCoefficients(primarySpecies_,
                                                 aqComplexRxns_);
  for (std::vector<Species>::iterator i = primarySpecies_.begin();
       i != primarySpecies_.end(); i++)
    i->update();

}

void Beaker::initializeMolalities(double initial_molality) 
{
  for (std::vector<Species>::iterator i=primarySpecies_.begin();
       i != primarySpecies_.end(); i++)
    i->molality(initial_molality);
} // end initializeMolalities

void Beaker::initializeMolalities(std::vector<double> initial_molalities) 
{
  if (initial_molalities.size() != primarySpecies_.size()) {
    std::cout << "Mismatch in size of initial_molalities array "
              << "and number of primarySpecies" << std::endl;
    exit(EXIT_SUCCESS);
  }

  // iterator doesnt seem to work then passing a vector entry - geh
  for (int i = 0; i < (int)primarySpecies_.size(); i++)
    primarySpecies_[i].molality(initial_molalities[i]);
} // end initializeMolalities()

void Beaker::updateEquilibriumChemistry(void)
{

  //    calculateActivityCoefficients(-1);

  // update primary species activities
  for (std::vector<Species>::iterator i = primarySpecies_.begin();
       i != primarySpecies_.end(); i++) {
    i->update();
  }
  // calculated seconday aqueous complex concentrations
  for (std::vector<AqueousEquilibriumComplex>::iterator i = aqComplexRxns_.begin();
       i != aqComplexRxns_.end(); i++) {
    i->update_kludge(primarySpecies_);
  }
  // calculate total component concentrations
  calculateTotal();

  // add equilibrium surface complexation here

} // end updateEquilibriumChemistry()

void Beaker::calculateTotal(std::vector<double> &total) 
{
  // add in primaries
  for (int i = 0; i < (int)total.size(); i++)
    total[i] = primarySpecies_[i].molality();

  // add in aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::iterator i = aqComplexRxns_.begin();
       i != aqComplexRxns_.end(); i++) {
    i->addContributionToTotal(total);
  }
  
  // scale by water density to convert to molarity
  for (int i = 0; i < (int)total.size(); i++)
    total[i] *= water_density_kg_L();
} // end calculateTotal()

void Beaker::calculateTotal(void) 
{
  calculateTotal(total_);
} // end calculateTotal()

void Beaker::calculateDTotal(Block *dtotal) 
{
  dtotal->zero();
  // derivative with respect to free-ion is 1.
  dtotal->setDiagonal(1.);

  // add in derviative of complex contribution with respect to free-ion
  for (std::vector<AqueousEquilibriumComplex>::iterator i = aqComplexRxns_.begin();
       i != aqComplexRxns_.end(); i++)
    i->addContributionToDTotal(primarySpecies_,dtotal);
  
  // scale by density of water
  dtotal->scale(water_density_kg_L()); 
} // end calculateDTotal()

void Beaker::calculateDTotal(void) 
{
  calculateDTotal(dtotal_);
}

void Beaker::updateKineticChemistry(void)
{
  // loop over general kinetic reactions and update effective rates
  for (std::vector<GeneralRxn>::iterator i = generalKineticRxns_.begin();
       i != generalKineticRxns_.end(); i++) {
    i->update_rates(primarySpecies_);
  }
  
  // add mineral saturation and rate calculations here
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->Update(primarySpecies_);
  }
  // add multirate kinetic surface complexation reaction quotient calculations 
  // here

} // end updateKineticChemistry()

void Beaker::addKineticChemistryToResidual(std::vector<double> &residual) 
{
  // loop over general kinetic reactions and add rates
  for (std::vector<GeneralRxn>::iterator i = generalKineticRxns_.begin();
       i != generalKineticRxns_.end(); i++)
         i->addContributionToResidual(residual, por_sat_den_vol());

  // add mineral mineral contribution to residual here.  units = mol/sec.
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->AddContributionToResidual(por_sat_den_vol(), &residual);
  }

  // add multirate kinetic surface complexation contribution to residual here.

} // end addKineticChemistryToResidual()

void Beaker::addKineticChemistryToJacobian(Block *J) 
{
  // loop over general kinetic reactions and add rates
  for (std::vector<GeneralRxn>::iterator i = generalKineticRxns_.begin();
       i != generalKineticRxns_.end(); i++)
         i->addContributionToJacobian(J, primarySpecies_, por_sat_den_vol());

  // add mineral mineral contribution to Jacobian here.  units = kg water/sec.
  for (std::vector<KineticRate*>::iterator rate = mineral_rates_.begin();
       rate != mineral_rates_.end(); rate++) {
    (*rate)->AddContributionToJacobian(primarySpecies_, por_sat_den_vol(), J);
  }

  // add multirate kinetic surface complexation contribution to Jacobian here.

} // end addKineticChemistryToJacobian()



void Beaker::addAccumulation(std::vector<double> &residual)
{
  addAccumulation(total_,residual);
}

void Beaker::addAccumulation(std::vector<double> total,
                             std::vector<double> &residual)
{
  // accumulation_coef = porosity*saturation*volume*1000./dt
  // units = (mol solute/L water)*(m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  //         (m^3 bulk)*(1000L water/m^3 water)/(sec) = mol/sec
  // 1000.d0 converts vol from m^3 -> L
  // all residual entries should be in mol/sec
  for (int i = 0; i < (int)total.size(); i++)
    residual[i] += accumulation_coef()*total[i];

  // add accumulation term for equilibrium sorption (e.g. Kd, surface 
  // complexation) here

} // end calculateAccumulation()

void Beaker::addAccumulationDerivative(Block *J)
{
  addAccumulationDerivative(J, dtotal_);
} // end calculateAccumulationDerivative()

void Beaker::addAccumulationDerivative(Block *J,
                                       Block *dtotal)
{
  // accumulation_coef = porosity*saturation*volume*1000./dt
  // units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*(m^3 bulk)/(sec)
  //         *(kg water/L water)*(1000L water/m^3 water) = kg water/sec
  // all Jacobian entries should be in kg water/sec
  J->addValues(dtotal,accumulation_coef());

  // add accumulation derivative term for equilibrium sorption 
  // (e.g. Kd, surface complexation) here

} // end calculateAccumulationDerivative()

void Beaker::calculateFixedAccumulation(std::vector<double> total,
                                        std::vector<double> &fixed_accumulation)
{
  for (int i = 0; i < (int)total.size(); i++)
    fixed_accumulation[i] = 0.;
  addAccumulation(total,fixed_accumulation);
} // end calculateAccumulation()

void Beaker::calculateResidual(std::vector<double> &residual, 
                               std::vector<double> fixed_accumulation)
{
  // subtract fixed porition
  for (int i = 0; i < ncomp(); i++)
    residual[i] = -fixed_accumulation[i];

  // accumulation adds in equilibrium chemistry
  addAccumulation(residual);

  // kinetic reaction contribution to residual
  addKineticChemistryToResidual(residual);


} // end calculateResidual()

void Beaker::calculateJacobian(Block *J)
{
  // must calculate derivatives with 
  calculateDTotal();

  // zero Jacobian
  J->zero();
  // add in derivatives for equilibrium chemistry
  addAccumulationDerivative(J);
  // add in derivatives for kinetic chemistry
  addKineticChemistryToJacobian(J);
} // end calculateJacobian()

void Beaker::scaleRHSAndJacobian(double *rhs, Block *J) 
{

  for (int i = 0; i < J->getSize(); i++) {
    double max = J->getRowAbsMax(i);
    if (max > 1.) {
      double scale = 1./max;
      rhs[i] *= scale;
      J->scaleRow(i,scale);
    }
  }
} // end scaleRHSAndJacobian()

void Beaker::scaleRHSAndJacobian(std::vector<double> &rhs, Block *J) 
{

  for (int i = 0; i < J->getSize(); i++) {
    double max = J->getRowAbsMax(i);
    if (max > 1.) {
      double scale = 1./max;
      rhs[i] *= scale;
      J->scaleRow(i,scale);
    }
  }
} // end scaleRHSAndJacobian()

void Beaker::updateMolalitiesWithTruncation(std::vector<double> &update, 
                                            std::vector<double> &prev_solution,
                                            double max_change) 
{
  // truncate the update to max_change
  for (int i = 0; i < ncomp(); i++) {
    if (update[i] > max_change) {
      update[i] = max_change;
    }
    else if (update[i] < -max_change) {
      update[i] = -max_change;
    }
    // store the previous solution
    prev_solution[i] = primarySpecies_[i].molality();
    // update primary species molalities (log formulation)
    double molality = prev_solution[i] * std::exp(-update[i]);
    primarySpecies_[i].molality(molality);
  }
} // end updateMolalitiesWithTruncation()

double Beaker::calculateMaxRelChangeInMolality(std::vector<double> prev_molal)
{
  double max_rel_change = 0.0;
  for (int i = 0; i < ncomp(); i++) {
    double delta = fabs(primarySpecies_[i].molality()-prev_molal[i]) / prev_molal[i];
    max_rel_change = delta > max_rel_change ? delta : max_rel_change;
  }
  return max_rel_change;
} // end calculateMaxRelChangeInMolality()

void Beaker::solveLinearSystem(Block *A, std::vector<double> &b) {
      // LU direct solve
    double D;
    // allocate pivoting array for LU
    int *indices = new int[ncomp()];
    ludcmp(A->getValues(),ncomp(),indices,&D);
    lubksb(A->getValues(),ncomp(),indices,b);
} // end solveLinearSystem

int Beaker::ReactionStep(std::vector<double> &total, 
			 const Beaker::BeakerParameters& parameters,
			 double dt)
{

  // update class paramters
  // water_density [kg/m^3]
  // volume [m^3]
  updateParameters(parameters, dt);

  // store current molalities
  for (int i = 0; i < ncomp(); i++)
    prev_molal[i] = primarySpecies_[i].molality();

  // initialize to a large number (not necessary, but safe)
  double max_rel_change = 1.e20;
  // iteration counter
  unsigned int num_iterations = 0;

  // lagging activity coefficients by a time step in this case
  updateActivityCoefficients();

  // calculate portion of residual at time level t
  calculateFixedAccumulation(total,fixed_accumulation);
                        
  do {

    // update equilibrium and kinetic chemistry (rates, ion activity 
    // products, etc.)
    updateEquilibriumChemistry();
    updateKineticChemistry();

    // units of residual: mol/sec
    calculateResidual(residual,fixed_accumulation);
    // units of Jacobian: kg water/sec
    calculateJacobian(J);
    // therefore, units of solution: mol/kg water (change in molality)

    for (int i = 0; i<ncomp(); i++)
      rhs[i] = residual[i];

    if (verbosity() == kDebugBeaker) {
      print_linear_system("before scale",J,rhs);
    }

    // scale the Jacobian
    scaleRHSAndJacobian(rhs,J);

    if (verbosity() == kDebugBeaker) {
      print_linear_system("after scale",J,rhs);
    }
    // for derivatives with respect to ln concentration, scale columns
    // by primary species concentrations
    for (int i = 0; i<ncomp(); i++)
      J->scaleColumn(i,primarySpecies_[i].molality());

    if (verbosity() == kDebugBeaker) {
      print_linear_system("before solve",J,rhs);
    }

    // call solver
    solveLinearSystem(J,rhs);

    if (verbosity() == kDebugBeaker) {
      print_linear_system("after solve",NULL,rhs);
    }

    // units of solution: mol/kg water (change in molality)
    // calculate update truncating at a maximum of 5 in nat log space
    updateMolalitiesWithTruncation(rhs,prev_molal,5.);
    // calculate maximum relative change in concentration over all species
    max_rel_change = calculateMaxRelChangeInMolality(prev_molal);

    if (verbosity() >= kDebugBeaker) {
      for (int i = 0; i < ncomp(); i++)
        std::cout << primarySpecies_[i].name() << " " << 
                  primarySpecies_[i].molality() << " " << total_[i] << "\n";
    }

    num_iterations++;

    // exit if maximum relative change is below tolerance
  } while (max_rel_change > tolerance() && num_iterations < max_iterations());

  // update total concentrations
  calculateTotal();
  for (int i = 0; i < ncomp(); i++)
    total[i] = total_[i];

  return num_iterations;

} // end ReactionStep()


// if no water density provided, default is 1000.0 kg/m^3
int Beaker::speciate(std::vector<double> target_total, const double water_density)
{
  double speciation_tolerance = 1.e-12;
  // water_density is in kg/m^3
  this->water_density_kg_m3(water_density);

  // initialize free-ion concentration s
  initializeMolalities(1.e-9);

  // store current molalities
  for (int i = 0; i < ncomp(); i++)
    prev_molal[i] = primarySpecies_[i].molality();

  double max_rel_change;
  unsigned int num_iterations = 0;
  bool calculate_activity_coefs = false;

  do {
    
    updateActivityCoefficients();
    updateEquilibriumChemistry();
    calculateDTotal();

    // calculate residual
    // units of residual: mol/sec
    for (int i = 0; i < ncomp(); i++)
      residual[i] = total_[i] - target_total[i];

    // add derivatives of total with respect to free to Jacobian
    // units of Jacobian: kg water/sec
    J->zero();
    calculateDTotal();
    J->addValues(0,0,dtotal_);

    // calculate residual
    for (int i = 0; i < ncomp(); i++)
      residual[i] = total_[i] - target_total[i];

    for (int i = 0; i < ncomp(); i++)
      rhs[i] = residual[i];

    if (verbosity() == kDebugBeaker) {
      print_linear_system("before scale",J,rhs);
    }

    // scale the Jacobian
    scaleRHSAndJacobian(rhs,J);

    if (verbosity() == kDebugBeaker) {
      print_linear_system("after scale",J,rhs);
    }
    // for derivatives with respect to ln concentration, scale columns
    // by primary species concentrations
    for (int i = 0; i<ncomp(); i++)
      J->scaleColumn(i,primarySpecies_[i].molality());

    if (verbosity() == kDebugBeaker) {
      print_linear_system("before solve",J,rhs);
    }

    // call solver
    solveLinearSystem(J,rhs);

    if (verbosity() == kDebugBeaker) {
      print_linear_system("after solve",NULL,rhs);
    }

    // calculate update truncating at a maximum of 5 in log space
    updateMolalitiesWithTruncation(rhs,prev_molal,5.);
    // calculate maximum relative change in concentration over all species
    max_rel_change = calculateMaxRelChangeInMolality(prev_molal);

    if (verbosity() == kDebugBeaker) {
      for (int i = 0; i < ncomp(); i++) {
      	std::cout << primarySpecies_[i].name() << " " 
		      << primarySpecies_[i].molality() << " " << total_[i] << "\n";
      }
    }

    num_iterations++;

    // if max_rel_change small enough, turn on activity coefficients
    if (max_rel_change < speciation_tolerance) calculate_activity_coefs = true;

    // exist if maximum relative change is below tolerance
  } while (max_rel_change > speciation_tolerance && 
	   num_iterations < max_iterations() && 
	   !calculate_activity_coefs);

  if (verbosity() > 1) {
    std::cout << "Beaker::speciate num_iterations :" << num_iterations << std::endl;
  }
  return num_iterations;
}  // end speciate()

void Beaker::display(void) const
{
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
} // end display()

void Beaker::print_results(void) const
{
  // output for testing purposes
  std::cout << std::endl;
  std::cout << "----- Solution ----------------------" << std::endl;
  std::cout << "Primary Species ---------------------\n";
  for (int i = 0; i < ncomp(); i++) {
    std::cout << "   " << primarySpecies_[i].name() << std::endl;
    std::cout << "        Total: " << total_[i] << std::endl;
    std::cout << "     Free-Ion: " << primarySpecies_[i].molality() << std::endl;
    std::cout << "Activity Coef: " << primarySpecies_[i].act_coef() << std::endl;
    std::cout << "     Activity: " << primarySpecies_[i].activity() << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Secondary Species -------------------\n";
  for (int i = 0; i < (int)aqComplexRxns_.size(); i++) {
    std::cout << "   " << aqComplexRxns_[i].name() << std::endl;
    std::cout << "     Free-Ion: " << aqComplexRxns_[i].molality() << std::endl;
    std::cout << "Activity Coef: " << aqComplexRxns_[i].act_coef() << std::endl;
    std::cout << "     Activity: " << aqComplexRxns_[i].activity() << std::endl;
  }
  std::cout << "-------------------------------------\n";
  std::cout << std::endl;

} // end print_results()

void Beaker::print_results(double time) const
{
  if (time < 1.e-40) {
    std::cout << "Time\t";
    for (int i = 0; i < ncomp(); i++) {
      std::cout << primarySpecies_[i].name() << " (total)\t";
      std::cout << primarySpecies_[i].name() << " (free-ion)\t";
    }
    std::cout << std::endl;
  }
  // output for testing purposes
  std::cout << time << "\t";
  for (int i = 0; i < ncomp(); i++) {
    std::cout << total_[i] << "\t";
    std::cout << primarySpecies_[i].molality() << "\t";
  }
  std::cout << std::endl;
} // end print_results()

void Beaker::print_linear_system(string s, Block *A, 
                                 std::vector<double> vector) {
  std::cout << s << endl;
  for (int i = 0; i < (int)vector.size(); i++)
    std::cout << "RHS: " << primarySpecies_[i].name() << " " << vector[i] << std::endl;
  if (A) A->print();
} // end print_linear_system()


