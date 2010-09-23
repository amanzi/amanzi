#include "Beaker.hpp"

Beaker::Beaker() 
{
  // end Beaker() constructor
}

Beaker::~Beaker() 
{
  // end Beaker destructor
}

void Beaker::setup(std::vector<double> *total) 
{
  
  /* end setup() */
}

void Beaker::addPrimarySpecies(Species s) 
{
  primarySpecies_.push_back(s);
}

void Beaker::addAqueousEquilibriumComplex(AqueousEquilibriumComplex c) 
{
  aqComplexRxns_.push_back(c);
}

void Beaker::updateParameters(double por, double sat, double vol, double dtt)
{
  porosity(por);
  saturation(sat);
  volume(vol);
  dt(dtt);
  psv1k_t(por,sat,vol,dtt);
}

void Beaker::initializeMolalities(double initial_molality) 
{
  for (std::vector<Species>::iterator i=primarySpecies_.begin();
       i!=primarySpecies_.end(); i++)
    i->set_molality(initial_molality);
}

void Beaker::initializeMolalities(std::vector<double> initial_molalities) 
{
  if (initial_molalities.size() != primarySpecies_.size()) {
    std::cout << "Mismatch in size of initial_molalities array "
              << "and number of primarySpecies" << std::endl;
    exit(EXIT_SUCCESS);
  }

//  for (std::vector<Species>::iterator i=primarySpecies_.begin();
//       i!=primarySpecies_.end(); i++)
//    i->set_molality(initial_molalities[i]);
  for (int i=0; i<(int)primarySpecies_.size(); i++)
    primarySpecies_[i].set_molality(initial_molalities[i]);
}

void Beaker::updateChemistry(void)
{
  for (std::vector<Species>::iterator i = primarySpecies_.begin();
       i != primarySpecies_.end(); i++) {
    i->update();
  }
  for (std::vector<AqueousEquilibriumComplex>::iterator i = aqComplexRxns_.begin();
       i != aqComplexRxns_.end(); i++) {
    i->update(primarySpecies_);
  }
}

void Beaker::calculateTotal(std::vector<double> &total) 
{
  // add in primaries
  for (int i = 0; i < (int)total.size(); i++)
    total[i] = primarySpecies_[i].get_molality();

  // add in aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::iterator i = aqComplexRxns_.begin();
       i != aqComplexRxns_.end(); i++) {
    i->addContributionToTotal(total);
  }

}

void Beaker::calculateDTotal(Block *dtotal) 
{

  dtotal->zero();
  // derivative with respect to free-ion is 1.
  dtotal->setDiagonal(1.);

  // add in derviative of complex contribution with respect to free-ion
  for (std::vector<AqueousEquilibriumComplex>::iterator i=aqComplexRxns_.begin();
       i!=aqComplexRxns_.end(); i++)
    i->addContributionToDTotal(primarySpecies_,dtotal);
//dtotal->scale(den_kg_per_L); scale by density of water

}

void Beaker::scaleRHSAndJacobian(double *rhs, Block *J) 
{

  for (int i=0; i<J->getSize(); i++) {
    double max = J->getRowAbsMax(i);
    if (max > 1.) {
      double scale = 1./max;
      rhs[i] *= scale;
      J->scaleRow(i,scale);
    }
  }

}

void Beaker::scaleRHSAndJacobian(std::vector<double> &rhs, Block *J) 
{

  for (int i=0; i<J->getSize(); i++) {
    double max = J->getRowAbsMax(i);
    if (max > 1.) {
      double scale = 1./max;
      rhs[i] *= scale;
      J->scaleRow(i,scale);
    }
  }

}

void Beaker::calculateAccumulation(std::vector<double> total,
                                   std::vector<double> &residual)
{
  // psv1k_t = porosity*saturation*volume*1000./dt
  // units = (mol solute/L water)*(m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  //         (m^3 bulk)*(1000L water/m^3 water)/(sec) = mol/sec
  // 1000.d0 converts vol from m^3 -> L
  // all residual entries should be in mol/sec
  for (int i=0; i<(int)total.size(); i++)
    residual[i] = psv1k_t()*total[i];
}

void Beaker::calculateAccumulationDerivative(Block *dtotal,
                                             Block *J)
{
  // psv1k_t = porosity*saturation*volume*1000./dt
  // units = (m^3 por/m^3 bulk)*(m^3 water/m^3 por)*(m^3 bulk)/(sec)
  //         *(kg water/L water)*(1000L water/m^3 water) = kg water/sec
  // all Jacobian entries should be in kg water/sec
  // note that setValues() overwrites all values...no need to zero
  J->setValues(dtotal,psv1k_t());
}


void Beaker::updateMolalitiesWithTruncation(std::vector<double> &update, 
                                            std::vector<double> &prev_solution,
                                            double max_change) 
{
  for (int i = 0; i < ncomp(); i++) {
    if (update[i] > max_change) update[i] = max_change;
    else if (update[i] < -max_change) update[i] = -max_change;
    prev_solution[i] = primarySpecies_[i].get_molality();
    primarySpecies_[i].set_molality(prev_solution[i]*exp(-update[i]));
  }
}

void Beaker::calculateMaxRelChangeInMolality(std::vector<double> prev_molal, 
                                             double &max_rel_change)
{
  max_rel_change = 0.;
  for (int i=0; i<ncomp(); i++) {
    double delta = fabs(primarySpecies_[i].get_molality()-prev_molal[i])/prev_molal[i];
    max_rel_change = delta > max_rel_change ? delta : max_rel_change;
  }
}

void Beaker::solveLinearSystem(Block *A, std::vector<double> &b) {
      // LU direct solve
    double D;
    // allocate pivoting array for LU
    int *indices = new int[ncomp()];
    ludcmp(A->getValues(),ncomp(),indices,&D);
    lubksb(A->getValues(),ncomp(),indices,b);
}

int Beaker::react(std::vector<double> &total, double volume, 
                  double porosity, double saturation, double dt)
{

  // update class paramters
  updateParameters(porosity,saturation,volume,dt);
  initializeMolalities(1.e-9);

  double tolerance = 1.e-12;

  Block *dtotal = new Block(ncomp());
  Block *J = new Block(ncomp());

  std::vector<double> fixed_residual(ncomp());
  std::vector<double> residual(ncomp());
  std::vector<double> rhs(ncomp());
  std::vector<double> prev_molal(ncomp());

  double max_rel_change;
  int num_iterations = 0;

  // calculate portion of residual at time level t
  calculateAccumulation(total,fixed_residual);
                        
  do {

//    calculateActivityCoefficients(-1);
    updateChemistry();
    calculateTotal(total);
    calculateDTotal(dtotal);

    calculateAccumulation(total,residual);
    // add derivatives of total with respect to free to Jacobian
    // calculateAccumulationDerivative() overwrites the entire Jacobian
    // i.e. no need to zero
    calculateAccumulationDerivative(dtotal,J);
  
    // subtract fixed porition
    for (int i=0; i<ncomp(); i++)
      residual[i] -= fixed_residual[i];

    // add additional reactions here
    // need contributions to residual and Jacobian
    
    // equilibrium surface complexation placeholder
    // calculateAccumulationSorb()
    // calculateAccumulationDerivSorb()

    if (verbose() == 3) {
      cout << "before scale\n";
      for (int i=0; i<ncomp(); i++)
          std::cout << "Res: " << primarySpecies_[i].get_name() << " " << 
                    residual[i] << "\n";
        J->print();
    }
    // scale the Jacobian
    for (int i=0; i<ncomp(); i++)
      rhs[i] = residual[i];
    scaleRHSAndJacobian(rhs,J);

    if (verbose() == 3) {
      cout << "after scale\n";
      for (int i=0; i<ncomp(); i++)
          std::cout << "RHS: " << primarySpecies_[i].get_name() << " " << 
                    rhs[i] << "\n";
      J->print();
    }
    // for derivatives with respect to ln concentration, scale columns
    // by primary species concentrations
    for (int i=0; i<ncomp(); i++)
      J->scaleColumn(i,primarySpecies_[i].get_molality());

    if (verbose() == 3) {
      cout << "before solve\n";
      for (int i=0; i<ncomp(); i++)
        std::cout << "RHS: " << primarySpecies_[i].get_name() << " " << 
                  rhs[i] << "\n";
      J->print();
    }

    // call solver
    solveLinearSystem(J,rhs);
    if (verbose() >= 3) {
      for (int i=0; i<ncomp(); i++)
        std::cout << "Update: " << primarySpecies_[i].get_name() << " " << 
                  rhs[i] << "\n";
    }

    // calculate update truncating at a maximum of 5 in log space
    updateMolalitiesWithTruncation(rhs,prev_molal,5.);
    // calculate maximum relative change in concentration over all species
    calculateMaxRelChangeInMolality(prev_molal,max_rel_change);

    if (verbose() >= 2) {
      for (int i=0; i<ncomp(); i++)
        std::cout << primarySpecies_[i].get_name() << " " << 
                  primarySpecies_[i].get_molality() << " " << total[i] << "\n";
    }

    num_iterations++;

    // exit if maximum relative change is below tolerance
  } while (max_rel_change > tolerance);

  totals_.resize(ncomp());
  for (int i = 0; i < ncomp(); i++) {
    totals_[i] = total[i];
  }

  // free up memory
  delete J;
  delete dtotal;

  return num_iterations;

}

int Beaker::speciate(std::vector<double> target_totals)
{

  double speciation_tolerance = 1.e-12;

  // initialize free-ion concentration s
  initializeMolalities(1.e-9);


  // allocate arrays for Newton-Raphson
  std::vector<double> totals(ncomp(), 0.0);
  Block *dtotal = new Block(ncomp());
  Block *J = new Block(ncomp());
  double *residual = new double[ncomp()];
  double *rhs = new double[ncomp()];
  double *prev_molal = new double[ncomp()];
  double *update = new double[ncomp()];

  // allocate pivoting array for LU
  int *indices = new int[ncomp()];

  double max_rel_change;
  int num_iterations = 0;

  do {
    
    //    calculateActivityCoefficients(-1);
    updateChemistry();
    calculateTotal(totals);
    calculateDTotal(dtotal);

    // add derivatives of total with respect to free to Jacobian
    J->zero();
    J->addValues(0,0,dtotal);

    // calculate residual
    for (int i = 0; i < ncomp(); i++) {
      residual[i] = totals[i] - target_totals[i];
    }

    if (verbose() == 3) {
      std::cout << "before scale\n";
      J->print();
    }

    // scale the Jacobian
    for (int i = 0; i < ncomp(); i++) {
      rhs[i] = residual[i];
    }
    scaleRHSAndJacobian(rhs, J);

    if (verbose() == 3) {
      std::cout << "after scale\n";
      J->print();
    }

    // for derivatives with respect to ln concentration, scale columns
    // by primary species concentrations
    for (int i = 0; i < ncomp(); i++) {
      J->scaleColumn(i, primarySpecies_[i].get_molality());
    }

    if (verbose() == 3) {
      std::cout << "before solve\n";
      J->print();
    }

    // LU direct solve
    double D;
    ludcmp(J->getValues(), ncomp(), indices, &D);
    lubksb(J->getValues(), ncomp(), indices, rhs);

    // the following two sections still need to be encapsulated in
    // function calls.

    // calculate update truncating at a maximum of 5 in log space
    for (int i = 0; i < ncomp(); i++) {
      update[i] = rhs[i] > 0. ? 
        (rhs[i] > 5. ? 5. : rhs[i]) : (rhs[i] < -5. ? -5. : rhs[i]);
      prev_molal[i] = primarySpecies_[i].get_molality();
      primarySpecies_[i].set_molality(prev_molal[i]*exp(-update[i]));
    }

    // calculate maximum relative change in concentration over all species
    max_rel_change = 0.;
    for (int i = 0; i < ncomp(); i++) {
      double delta = fabs(primarySpecies_[i].get_molality() - prev_molal[i]) / prev_molal[i];
      max_rel_change = delta > max_rel_change ? delta : max_rel_change;
    }

    if (verbose() == 3) {
      for (int i = 0; i < ncomp(); i++) {
      	std::cout << primarySpecies_[i].get_name() << " " 
		      << primarySpecies_[i].get_molality() << " " << totals[i] << "\n";
      }
    }

    num_iterations++;

    // exist if maximum relative change is below tolerance
  } while (max_rel_change > speciation_tolerance);

  // free up memory
  delete J;
  delete dtotal;
  delete [] residual;
  delete [] rhs;
  delete [] update;
  delete [] prev_molal;
  delete [] indices;

  totals_.resize(ncomp());
  for (int i = 0; i < ncomp(); i++) {
    totals_[i] = totals[i];
  }

  if (verbose() > 1) {
    std::cout << "Beaker::speciate num_iterations :" << num_iterations << std::endl;
  }
  return num_iterations;
  // end speciate
}

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
  // end display()
}

void Beaker::print_results(void) const
{
  // output for testing purposes
  std::cout << std::endl;
  std::cout << "----- Solution ----------------------" << std::endl;
  std::cout << "Primary Species ---------------------\n";
  for (int i = 0; i < ncomp(); i++) {
    std::cout << "  " << primarySpecies_[i].get_name() << std::endl;
    std::cout << "       Total: " << totals_[i] << std::endl;
    std::cout << "    Free-Ion: " << primarySpecies_[i].get_molality() << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Secondary Species -------------------\n";
  for (int i = 0; i < (int)aqComplexRxns_.size(); i++) {
    std::cout << "  " << aqComplexRxns_[i].get_name() << std::endl;
    std::cout << "    Free-Ion: " << aqComplexRxns_[i].get_molality() << std::endl;
  }
  std::cout << "-------------------------------------\n";
  std::cout << std::endl;

  // end print_results
}

