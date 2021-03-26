/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre
*/

#include "activity_model_unit.hh"

#include <cmath>

#include <iostream>

namespace Amanzi {
namespace AmanziChemistry {

double ActivityModelUnit::Evaluate(const Species& species) {
  static_cast<void>(species);
  // log(gamma_i) = 0.0, gamma_i = 1.0
  return 1.0;
}


void ActivityModelUnit::EvaluateVector(
    const std::vector<Species>& prim, 
    const std::vector<AqueousEquilibriumComplex>& sec,
    std::vector<double>* gamma, 
    double* actw) {
  // double r1(1.0);
  for (std::vector<double>::iterator i = gamma->begin(); i != gamma->end(); ++i) {
    (*i) = 1.0;
  }
  *actw = 1.0;
}


void ActivityModelUnit::Display(void) const {
  std::cout << "Activity Model: unit activity coefficients (gamma = 1.0)." << std::endl;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi

