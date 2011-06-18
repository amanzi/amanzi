/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "activity_model_unit.hh"

#include <cmath>

#include <iostream>

namespace amanzi {
namespace chemistry {

ActivityModelUnit::ActivityModelUnit()
    : ActivityModel() {
}  // end ActivityModelUnit constructor


ActivityModelUnit::~ActivityModelUnit() {
}  // end ActivityModelUnit destructor

double ActivityModelUnit::Evaluate(const Species& species) {
  static_cast<void>(species);
  // log(gamma_i) = 0.0, gamma_i = 1.0

  return 1.0;
}  // end Evaluate()

void ActivityModelUnit::Display(void) const {
  std::cout << "Activity Model: unit activity coefficients (gamma = 1.0)."
            << std::endl;
}  // end Display()

}  // namespace chemistry
}  // namespace amanzi
