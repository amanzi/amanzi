/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cmath>

#include <iostream>

#include "ActivityModelUnit.hpp"

ActivityModelUnit::ActivityModelUnit()
  : ActivityModel()
{
}  // end ActivityModelUnit constructor


ActivityModelUnit::~ActivityModelUnit()
{
}  // end ActivityModelUnit destructor

double ActivityModelUnit::evaluate(const Species& species)
{
  static_cast<void>(species);
  // log(gamma_i) = 0.0, gamma_i = 1.0

  return 1.0;
}  // end calculateActivityCoefficient()

void ActivityModelUnit::display(void) const
{
  std::cout << "Using unit activity coefficients (gamma = 1.0)." << std::endl;
}  // end display()
