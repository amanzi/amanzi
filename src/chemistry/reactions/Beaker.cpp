#include <cstdlib>

#include "Beaker.hpp"

Beaker::Beaker() 
{
  dtotal_ = NULL;
  J = NULL;
} // end Beaker() constructor

Beaker::~Beaker() 
{
  if (dtotal_) delete dtotal_;
  dtotal_ = NULL;
  if (J) delete J;
  J = NULL;
} // end Beaker destructor

void Beaker::resize() {
  if (dtotal_) delete dtotal_;
  dtotal_ = new Block(ncomp());
  fixed_residual.resize(ncomp());
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

void Beaker::setup(std::vector<double> &total) 
{
} // end setup()

void Beaker::addPrimarySpecies(Species s) 
{
  primarySpecies_.push_back(s);
} // end addPrimarySpecies()

void Beaker::addAqueousEquilibriumComplex(AqueousEquilibriumComplex c) 
{
  aqComplexRxns_.push_back(c);
} // end addAqueousEquilibriumComplex()

void Beaker::updateParameters(double por, double sat, double vol, double dtt)
{
  porosity(por);
  saturation(sat);
  volume(vol);
  dt(dtt);
  accumulation_coef(por,sat,vol,dtt);
} // end updateParameters()

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
  for (std::vector<Species>::iterator i = primarySpecies_.begin();
       i != primarySpecies_.end(); i++) {
    i->update();
  }
  for (std::vector<AqueousEquilibriumComplex>::iterator i = aqComplexRxns_.begin();
       i != aqComplexRxns_.end(); i++) {
    i->update(primarySpecies_);
  }
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
} // end calculateTotal()

void Beaker::calculateTotal() 
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
//dtotal->scale(den_kg_per_L); scale by density of water
} // end calculateDTotal()

void Beaker::calculateDTotal() 
{
  calculateDTotal(dtotal_);
}

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

void Beaker::calculateAccumulation(std::vector<double> &residual)
{
  calculateAccumulation(total_,residual);
}

void Beaker::calculateAccumulation(std::vector<double> total,
                                   std::vector<double> &residual)
{
  // psv1k_t = porosity*saturation*volume*1000./dt
  // units = (mol solute/L water)*(m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  //         (m^3 bulk)*(1000L water/m^3 water)/(sec) = mol/sec
  // 1000.d0 converts vol from m^3 -> L
  // all residual entries should be in mol/sec
  for (int i = 0; i < (int)total.size(); i++)
    residual[i] = accumulation_coef()*total[i];
} // end calculateAccumulation()

void Beaker::calculateAccumulationDerivative(Block *J)
{
  calculateAccumulationDerivative(dtotal_,J);
} // end calculateAccumulationDerivative()

void Beaker::calculateAccumulationDerivative(Block *dtotal,
                                             Block *J)
{
  // psv1k_t = porosity*saturation*volume*1000./dt
  // units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*(m^3 bulk)/(sec)
  //         *(kg water/L water)*(1000L water/m^3 water) = kg water/sec
  // all Jacobian entries should be in kg water/sec
  // note that setValues() overwrites all values...no need to zero
  J->setValues(dtotal,accumulation_coef());
} // end calculateAccumulationDerivative()


void Beaker::updateMolalitiesWithTruncation(std::vector<double> &update, 
                                            std::vector<double> &prev_solution,
                                            double max_change) 
{
  for (int i = 0; i < ncomp(); i++) {
    if (update[i] > max_change) {
      update[i] = max_change;
    }
    else if (update[i] < -max_change) {
      update[i] = -max_change;
    }
    prev_solution[i] = primarySpecies_[i].molality();
    primarySpecies_[i].molality(prev_solution[i]*exp(-update[i]));
  }
} // end updateMolalitiesWithTruncation()

void Beaker::calculateMaxRelChangeInMolality(std::vector<double> prev_molal, 
                                             double &max_rel_change)
{
  max_rel_change = 0.;
  for (int i = 0; i < ncomp(); i++) {
    double delta = fabs(primarySpecies_[i].molality()-prev_molal[i]) / prev_molal[i];
    max_rel_change = delta > max_rel_change ? delta : max_rel_change;
  }
} // end calculateMaxRelChangeInMolality()

void Beaker::solveLinearSystem(Block *A, std::vector<double> &b) {
      // LU direct solve
    double D;
    // allocate pivoting array for LU
    int *indices = new int[ncomp()];
    ludcmp(A->getValues(),ncomp(),indices,&D);
    lubksb(A->getValues(),ncomp(),indices,b);
} // end solveLinearSystem

int Beaker::react(std::vector<double> &total, double volume, 
                  double porosity, double saturation, double dt)
{

  // update class paramters
  updateParameters(porosity,saturation,volume,dt);
  initializeMolalities(1.e-9);

  double tolerance = 1.e-12;

  double max_rel_change;
  int num_iterations = 0;

  // calculate portion of residual at time level t
  calculateAccumulation(total,fixed_residual);
                        
  do {

//    calculateActivityCoefficients(-1);
    updateEquilibriumChemistry();
    calculateTotal();
    calculateDTotal();

    calculateAccumulation(residual);
    // add derivatives of total with respect to free to Jacobian
    // calculateAccumulationDerivative() overwrites the entire Jacobian
    // i.e. no need to zero
    calculateAccumulationDerivative(J);
  
    // subtract fixed porition
    for (int i = 0; i < ncomp(); i++)
      residual[i] -= fixed_residual[i];

    // add additional reactions here
    // need contributions to residual and Jacobian
    
    // equilibrium surface complexation placeholder
    // calculateAccumulationSorb()
    // calculateAccumulationDerivSorb()

    if (verbose() == 3) {
      print_linear_system("before scale",J,rhs);
    }

    // scale the Jacobian
    for (int i = 0; i<ncomp(); i++)
      rhs[i] = residual[i];
    scaleRHSAndJacobian(rhs,J);

    if (verbose() == 3) {
      print_linear_system("after scale",J,rhs);
    }
    // for derivatives with respect to ln concentration, scale columns
    // by primary species concentrations
    for (int i = 0; i<ncomp(); i++)
      J->scaleColumn(i,primarySpecies_[i].molality());

    if (verbose() == 3) {
      print_linear_system("before solve",J,rhs);
    }

    // call solver
    solveLinearSystem(J,rhs);

    if (verbose() >= 3) {
      print_linear_system("after solve",NULL,rhs);
    }

    // calculate update truncating at a maximum of 5 in log space
    updateMolalitiesWithTruncation(rhs,prev_molal,5.);
    // calculate maximum relative change in concentration over all species
    calculateMaxRelChangeInMolality(prev_molal,max_rel_change);

    if (verbose() >= 2) {
      for (int i = 0; i < ncomp(); i++)
        std::cout << primarySpecies_[i].name() << " " << 
                  primarySpecies_[i].molality() << " " << total[i] << "\n";
    }

    num_iterations++;

    // exit if maximum relative change is below tolerance
  } while (max_rel_change > tolerance);

  total_.resize(ncomp());
  for (int i = 0; i < ncomp(); i++) {
    total_[i] = total[i];
  }

  return num_iterations;

} // end react()

int Beaker::speciate(std::vector<double> target_total)
{

  double speciation_tolerance = 1.e-12;

  // initialize free-ion concentration s
  initializeMolalities(1.e-9);

  double max_rel_change;
  int num_iterations = 0;

  do {
    
    //    calculateActivityCoefficients(-1);
    updateEquilibriumChemistry();
    calculateTotal();
    calculateDTotal();

    // add derivatives of total with respect to free to Jacobian
    J->zero();
    J->addValues(0,0,dtotal_);

    // calculate residual
    for (int i = 0; i < ncomp(); i++)
      residual[i] = total_[i] - target_total[i];

    if (verbose() == 3) {
      print_linear_system("before scale",J,rhs);
    }

    // scale the Jacobian
    for (int i = 0; i < ncomp(); i++)
      rhs[i] = residual[i];
    scaleRHSAndJacobian(rhs,J);

    if (verbose() == 3) {
      print_linear_system("after scale",J,rhs);
    }
    // for derivatives with respect to ln concentration, scale columns
    // by primary species concentrations
    for (int i = 0; i<ncomp(); i++)
      J->scaleColumn(i,primarySpecies_[i].molality());

    if (verbose() == 3) {
      print_linear_system("before solve",J,rhs);
    }

    // call solver
    solveLinearSystem(J,rhs);

    if (verbose() >= 3) {
      print_linear_system("after solve",NULL,rhs);
    }

    // calculate update truncating at a maximum of 5 in log space
    updateMolalitiesWithTruncation(rhs,prev_molal,5.);
    // calculate maximum relative change in concentration over all species
    calculateMaxRelChangeInMolality(prev_molal,max_rel_change);

    if (verbose() == 3) {
      for (int i = 0; i < ncomp(); i++) {
      	std::cout << primarySpecies_[i].name() << " " 
		      << primarySpecies_[i].molality() << " " << total_[i] << "\n";
      }
    }

    num_iterations++;

    // exist if maximum relative change is below tolerance
  } while (max_rel_change > speciation_tolerance);

  if (verbose() > 1) {
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
    std::cout << "  " << primarySpecies_[i].name() << std::endl;
    std::cout << "       Total: " << total_[i] << std::endl;
    std::cout << "    Free-Ion: " << primarySpecies_[i].molality() << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Secondary Species -------------------\n";
  for (int i = 0; i < (int)aqComplexRxns_.size(); i++) {
    std::cout << "  " << aqComplexRxns_[i].name() << std::endl;
    std::cout << "    Free-Ion: " << aqComplexRxns_[i].molality() << std::endl;
  }
  std::cout << "-------------------------------------\n";
  std::cout << std::endl;

} // end print_results()

void Beaker::print_linear_system(string s, Block *A, 
                                 std::vector<double> vector) {
  std::cout << s << endl;
  for (int i = 0; i < (int)vector.size(); i++)
    std::cout << primarySpecies_[i].name() << " " << vector[i] << std::endl;
  if (A) A->print();
} // end print_linear_system()


