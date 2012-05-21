/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "activity_model_debye_huckel.hh"
#include <cmath>
#include <iostream>

namespace amanzi {
namespace chemistry {

const double ActivityModelDebyeHuckel::debyeA = 0.5114;  // 25C
const double ActivityModelDebyeHuckel::debyeB = 0.3288;  // 25C
const double ActivityModelDebyeHuckel::debyeBdot = 0.0410;  // 25C

ActivityModelDebyeHuckel::ActivityModelDebyeHuckel()
    : ActivityModel() {
}  // end ActivityModelDebyeHuckel constructor


ActivityModelDebyeHuckel::~ActivityModelDebyeHuckel() {
}  // end ActivityModelDebyeHuckel destructor

double ActivityModelDebyeHuckel::Evaluate(const Species& species) {
  // log(gamma_i) = - A * z_i^2 * sqrt(I) / (1 + a0 * B * sqrt(I)) + Bdot * I
  double gamma(0.0);
  if (fabs(species.charge()) < 1.e-10) {
    // for now, neutral species activity = 1.
    gamma = 1.0;
  } else {
    double sqrt_I = std::sqrt(I_);

    double log_gamma = -debyeA * species.charge() * species.charge() * sqrt_I /
        (1.0 + species.ion_size_parameter() * debyeB * sqrt_I) +
        debyeBdot * I_;

    // bja: why not just std::pow(10.0, log_gamma)?
    gamma = std::exp(log_to_ln(log_gamma));
  }
  return gamma;
}  // end Evaluate()

void ActivityModelDebyeHuckel::EvaluateVector(
    const std::vector<Species>& prim,
    const std::vector<AqueousEquilibriumComplex>& sec,
    std::vector<double>* gamma, 
    double* actw) {
  int isp(-1);
  // For primary species
  for (std::vector<Species>::const_iterator i = prim.begin(); i != prim.end(); i++) {
    isp++;
    if (fabs((*i).charge()) < 1.0e-10) {
      gamma->at(isp) = 1.0;
    } else {
      gamma->at(isp) = Evaluate(*i);
    }
  }
  // For aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::const_iterator i = sec.begin(); i != sec.end(); i++) {
    isp++;
    if (fabs((*i).charge()) < 1.0e-10) {
      gamma->at(isp) = 1.0;
    } else {
      gamma->at(isp) = Evaluate(*i);
    }
  }
  *actw = 1.0;
}  // end evaluate

void ActivityModelDebyeHuckel::Display(void) const {
  std::cout << "Activity model: Debye-Huckel" << std::endl;
}  // end Display()

}  // namespace chemistry
}  // namespace amanzi
