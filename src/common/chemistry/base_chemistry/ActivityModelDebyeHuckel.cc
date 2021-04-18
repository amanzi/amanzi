/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre

  Class for activity calculations based on the Debye-Huckel B-dot equation.
   
  TODO(bandre): need to fix the name of this class to be
  DebyeHuckelBdot or something to distinguish it from a pure
  Debye-Huckel. Is it worth worrying about code reuse between
  debye-huckel and debye-huckel b-dot?
*/

#include <cmath>
#include <iostream>

#include "ActivityModelDebyeHuckel.hh"
#include "ChemistryUtilities.hh"

namespace Amanzi {
namespace AmanziChemistry {

namespace acu = Amanzi::AmanziChemistry::utilities;

const double ActivityModelDebyeHuckel::debyeA = 0.5114;  // 25C
const double ActivityModelDebyeHuckel::debyeB = 0.3288;  // 25C
const double ActivityModelDebyeHuckel::debyeBdot = 0.0410;  // 25C

/* ******************************************************************
* log(gamma_i) = -A z_i^2 sqrt(I) / (1 + a0 B sqrt(I)) + Bdot  I
****************************************************************** */
double ActivityModelDebyeHuckel::Evaluate(const Species& species)
{
  double gamma(0.0);
  if (std::fabs(species.charge()) < 1.0e-10) {
    // for now, neutral species activity = 1.
    gamma = 1.0;
  } else {
    double sqrt_I = std::sqrt(I_);

    double log_gamma = -debyeA * species.charge() * species.charge() * sqrt_I /
        (1.0 + species.ion_size_parameter() * debyeB * sqrt_I) +
        debyeBdot * I_;

    // KL: this is the only place where we can control diverging Newton method
    log_gamma = std::min(log_gamma, max_log_gamma_); 

    // bja: why not just std::pow(10.0, log_gamma)?
    gamma = std::exp(acu::log_to_ln(log_gamma));
  }
  return gamma;
}


void ActivityModelDebyeHuckel::EvaluateVector(
    const std::vector<Species>& primary_species,
    const std::vector<AqueousEquilibriumComplex>& secondary_species,
    std::vector<double>* gamma, 
    double* actw)
{
  int isp(0);

  for (auto it = primary_species.begin(); it != primary_species.end(); ++it, ++isp) {
    if (std::fabs((*it).charge()) < 1.0e-10) {
      gamma->at(isp) = 1.0;
    } else {
      gamma->at(isp) = Evaluate(*it);
    }
  }

  for (auto it = secondary_species.begin(); it != secondary_species.end(); ++it, ++isp) {
    if (std::fabs((*it).charge()) < 1.0e-10) {
      gamma->at(isp) = 1.0;
    } else {
      gamma->at(isp) = Evaluate(*it);
    }
  }

  *actw = 1.0;
}


void ActivityModelDebyeHuckel::Display() const {
  std::cout << "Activity model: Debye-Huckel" << std::endl;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
