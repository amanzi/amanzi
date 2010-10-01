/* -*-  mode: c++; c-default-style: "google-c-style"; indent-tabs-mode: nil -*- */
#include <cmath>

#include <iostream>

#include "ActivityModelDebyeHuckel.hpp"

const double ActivityModelDebyeHuckel::debyeA = 0.4939;
const double ActivityModelDebyeHuckel::debyeB = 0.3253;
const double ActivityModelDebyeHuckel::debyeBdot = 0.0374;

ActivityModelDebyeHuckel::ActivityModelDebyeHuckel()
    : ActivityModel()
{
}  // end ActivityModelDebyeHuckel constructor


ActivityModelDebyeHuckel::~ActivityModelDebyeHuckel()
{
}  // end ActivityModelDebyeHuckel destructor

double ActivityModelDebyeHuckel::evaluate(const Species& species)
{
  // log(gamma_i) = - A * z_i^2 * sqrt(I) / (1 + a0 * B * sqrt(I)) + Bdot * I

  double sqrt_I = std::sqrt(I_);

  double log_gamma = -debyeA * species.charge() * species.charge() * sqrt_I /
      (1.0 + species.ion_size_parameter() * debyeB * sqrt_I) +
      debyeBdot * I_;

  return std::exp(log_to_ln(log_gamma));
}  // end calculateActivityCoefficient()

void ActivityModelDebyeHuckel::display(void) const
{
  std::cout << "Using Debye-Huckel activity model." << std::endl;
}  // end display()
