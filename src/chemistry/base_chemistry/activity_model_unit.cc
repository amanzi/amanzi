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

void ActivityModelUnit::EvaluateVector (std::vector<double>& gamma, const std::vector<Species>* prim, const std::vector<AqueousEquilibriumComplex>* sec){
	const double r1(1.0e0);
	for (std::vector<double>::iterator i=gamma.begin(); i!=gamma.end(); i++) (*i)=r1;
} // end EvaluateVector


void ActivityModelUnit::Display(void) const {
  std::cout << "Activity Model: unit activity coefficients (gamma = 1.0)."
            << std::endl;
}  // end Display()

}  // namespace chemistry
}  // namespace amanzi
