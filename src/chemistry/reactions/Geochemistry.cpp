#include "Geochemistry.hpp"

Geochemistry::Geochemistry() {
}

void Geochemistry::setup(std::vector<double> *total) {
  
  /* end setup() */
}

void Geochemistry::addPrimarySpecies(Species s) {
  primarySpecies_.push_back(s);
}

void Geochemistry::addAqueousEquilibriumComplex(AqueousEquilibriumComplex c) {
  aqComplexRxns_.push_back(c);
}

void Geochemistry::initializeMolalities(double initial_molality) {
  for (std::vector<Species>::iterator i=primarySpecies_.begin();
       i!=primarySpecies_.end(); i++)
    i->set_molality(initial_molality);
}

void Geochemistry::updateChemistry(void) {
  for (std::vector<Species>::iterator i = primarySpecies_.begin();
       i != primarySpecies_.end(); i++) {
    i->update();
  }
  for (std::vector<AqueousEquilibriumComplex>::iterator i = aqComplexRxns_.begin();
       i != aqComplexRxns_.end(); i++) {
    i->update(primarySpecies_);
  }
}

void Geochemistry::calculateTotal(std::vector<double> &total) {

  // add in primaries
  for (int i = 0; i < total.size(); i++)
    total[i] = primarySpecies_[i].get_molality();

  // add in aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::iterator i = aqComplexRxns_.begin();
       i != aqComplexRxns_.end(); i++) {
    i->addContributionToTotal(total);
  }

}

void Geochemistry::calculateDTotal(Block *dtotal) {

  dtotal->zero();
  // derivative with respect to free-ion is 1.
  dtotal->setDiagonal(1.);

  // add in derviative of complex contribution with respect to free-ion
  for (std::vector<AqueousEquilibriumComplex>::iterator i=aqComplexRxns_.begin();
       i!=aqComplexRxns_.end(); i++)
    i->addContributionToDTotal(primarySpecies_,dtotal);
//dtotal->scale(den_kg_per_L); scale by density of water

}

void Geochemistry::scaleRHSAndJacobian(double *rhs, Block *J) {

  for (int i=0; i<J->getSize(); i++) {
    double max = J->getRowAbsMax(i);
    if (max > 1.) {
      double scale = 1./max;
      rhs[i] = rhs[i]*scale;
      J->scaleRow(i,scale);
    }
  }

}

int Geochemistry::speciate(std::vector<double> target_totals)
{

  double speciation_tolerance = 1.e-12;

  // initialize free-ion concentration s
  initializeMolalities(1.e-9);


  // allocate arrays for Newton-Raphson
  std::vector<double> totals(ncomp(), 0.0);
  Block *dtotal = new Block(ncomp());
  double *residual = new double[ncomp()];
  double *rhs = new double[ncomp()];
  double *prev_molal = new double[ncomp()];
  double *update = new double[ncomp()];
  Block *J = new Block(ncomp());

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
      cout << "before scale\n";
      J->print();
    }

    // scale the Jacobian
    for (int i = 0; i < ncomp(); i++) {
      rhs[i] = residual[i];
    }
    scaleRHSAndJacobian(rhs, J);

    if (verbose() == 3) {
      cout << "after scale\n";
      J->print();
    }

    // for derivatives with respect to ln concentration, scale columns
    // by primary species concentrations
    for (int i = 0; i < ncomp(); i++) {
      J->scaleColumn(i, primarySpecies_[i].get_molality());
    }

    if (verbose() == 3) {
      cout << "before solve\n";
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
    cout << "Geochemistry::speciate num_iterations :" << num_iterations << endl;
  }
  return num_iterations;
  // end speciate
}

void Geochemistry::print_results(void)
{
  // output for testing purposes
  cout << endl;
  cout << "Primary Species ---------------------\n";
  for (int i = 0; i < ncomp(); i++) {
    cout << "  " << primarySpecies_[i].get_name() << endl;
    cout << "       Total: " << totals_[i] << endl;
    cout << "    Free-Ion: " << primarySpecies_[i].get_molality() << endl;
  }
  cout << endl;
  cout << "Secondary Species -------------------\n";
  for (int i = 0; i < aqComplexRxns_.size(); i++) {
    cout << "  " << aqComplexRxns_[i].get_name() << endl;
    cout << "    Free-Ion: " << aqComplexRxns_[i].get_molality() << endl;
  }
  cout << "-------------------------------------\n";
  cout << endl;

  // end print_results
}

Geochemistry::~Geochemistry() {
}
