/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre
*/

#include <iostream>

#include "ActivityModelUnit.hh"

namespace Amanzi {
namespace AmanziChemistry {

double ActivityModelUnit::Evaluate(const Species& species) {
  return 1.0;
}


void ActivityModelUnit::EvaluateVector(
    const std::vector<Species>& primary_species, 
    const std::vector<AqueousEquilibriumComplex>& secondary_species,
    std::vector<double>* gamma, 
    double* actw)
{
  gamma->assign(gamma->size(), 1.0);
  *actw = 1.0;
}


void ActivityModelUnit::Display() const {
  std::cout << "Activity Model: unit activity coefficients (gamma = 1.0)." << std::endl;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi

