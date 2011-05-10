/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cmath>

#include <vector>
#include <iostream>

#include "ActivityModel.hpp"

ActivityModel::ActivityModel()
    : I_(0.0) {
}  // end ActivityModel constructor

ActivityModel::~ActivityModel() {
}  // end ActivityModel destructor

void ActivityModel::CalculateIonicStrength(
    const std::vector<Species>& primarySpecies,
    const std::vector<AqueousEquilibriumComplex>& secondarySpecies) {
  // I = 0.5 * sum_i(m_i*z_i^2)
  I_ = 0.;

  // primary species
  for (std::vector<Species>::const_iterator i = primarySpecies.begin();
       i != primarySpecies.end(); i++) {
    I_ += i->molality() * i->charge() * i->charge();
  }

  // secondary aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::const_iterator i = secondarySpecies.begin();
       i != secondarySpecies.end(); i++) {
    I_ += i->molality() * i->charge() * i->charge();
  }

  I_ *= 0.5;
}  // end CalculateIonicStrength()

void ActivityModel::CalculateActivityCoefficients(
    std::vector<Species>* primarySpecies,
    std::vector<AqueousEquilibriumComplex>* secondarySpecies) {
  // primary species
  for (std::vector<Species>::iterator i = primarySpecies->begin();
       i != primarySpecies->end(); i++) {
    double gamma = Evaluate(*i);
    i->act_coef(gamma);
  }

  // secondary aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::iterator i = secondarySpecies->begin();
       i != secondarySpecies->end(); i++) {
    double gamma = Evaluate(*i);
    i->act_coef(gamma);
  }
}  // end CalculateActivityCoefficients()

