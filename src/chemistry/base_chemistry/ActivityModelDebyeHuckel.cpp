/* -*-  mode: c++; c-default-style: "google-c-style"; indent-tabs-mode: nil -*- */
#include <cmath>

#include <iostream>

#include "ActivityModelDebyeHuckel.hpp"

const double ActivityModelDebyeHuckel::debyeA = 0.5114;
const double ActivityModelDebyeHuckel::debyeB = 0.3288;
const double ActivityModelDebyeHuckel::debyeBdot = 0.0410;

ActivityModelDebyeHuckel::ActivityModelDebyeHuckel()
    : ActivityModel()
{
}  // end ActivityModelDebyeHuckel constructor


ActivityModelDebyeHuckel::~ActivityModelDebyeHuckel()
{
}  // end ActivityModelDebyeHuckel destructor

double ActivityModelDebyeHuckel::Evaluate(const Species& species)
{
  // log(gamma_i) = - A * z_i^2 * sqrt(I) / (1 + a0 * B * sqrt(I)) + Bdot * I

  // for now, neutral species activity = 1.
  if (fabs(species.charge()) < 1.e-10) return 1.;

  double sqrt_I = std::sqrt(I_);

  double log_gamma = -debyeA * species.charge() * species.charge() * sqrt_I /
      (1.0 + species.ion_size_parameter() * debyeB * sqrt_I) +
      debyeBdot * I_;

  return std::exp(log_to_ln(log_gamma));
}  // end Evaluate()

void ActivityModelDebyeHuckel::Display(void) const
{
  std::cout << "Activity model: Debye-Huckel" << std::endl;
}  // end Display()
