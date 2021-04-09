/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre
*/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>

#include "activity_model.hh"

namespace Amanzi {
namespace AmanziChemistry {

ActivityModel::ActivityModel()
    : I_(0.0),
      Z_(0.0),
      M_(0.0),
      name_(""),
      num_species_(0),
      gamma_() {
}


void ActivityModel::Setup(
    const ActivityModelParameters& parameters,
    const std::vector<Species>& primary_species,
    const std::vector<AqueousEquilibriumComplex>& secondary_species) {
  ResizeGamma(primary_species.size() + secondary_species.size());
}


void ActivityModel::ResizeGamma(int size) {
  num_species_ = size;
  gamma_.resize(num_species_, 1.0);  
}


void ActivityModel::CalculateIonicStrength(
    const std::vector<Species>& primarySpecies,
    const std::vector<AqueousEquilibriumComplex>& secondarySpecies) {
  // I = 0.5 * sum_i(m_i*z_i^2)
  I_ = 0.0;

  // primary species
  for (auto it = primarySpecies.begin(); it != primarySpecies.end(); ++it) {
    I_ += it->molality() * it->charge() * it->charge();
  }

  // secondary aqueous complexes
  for (auto it = secondarySpecies.begin(); it != secondarySpecies.end(); ++it) {
    I_ += it->molality() * it->charge() * it->charge();
  }

  I_ *= 0.5;
}


void ActivityModel::CalculateSumAbsZ(
    const std::vector<Species>& primarySpecies,
    const std::vector<AqueousEquilibriumComplex>& secondarySpecies) {
  // Z = sum_i(m_i*abs(z_i))
  Z_ = 0.0;

  // primary species
  for (std::vector<Species>::const_iterator i = primarySpecies.begin();
       i != primarySpecies.end(); i++) {
    if (i->name() != "h2o" && i->name() != "H2O") {
      Z_ += i->molality() * std::abs(i->charge());
    }
  }

  // secondary aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::const_iterator i = secondarySpecies.begin();
       i != secondarySpecies.end(); i++) {
    if (i->name() != "h2o" && i->name() != "H2O") {
      Z_ += i->molality() * std::abs(i->charge());
    }
  }
}


void ActivityModel::CalculateSumC(
    const std::vector<Species>& primarySpecies,
    const std::vector<AqueousEquilibriumComplex>& secondarySpecies) {
  // Z = sum_i(m_i*abs(z_i))
  M_ = 0.0;

  // primary species
  for (std::vector<Species>::const_iterator i = primarySpecies.begin();
       i != primarySpecies.end(); i++) {
    if (i->name() != "h2o" && i->name() != "H2O") {
      M_ += i->molality();
    }
  }

  // secondary aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::const_iterator i = secondarySpecies.begin();
       i != secondarySpecies.end(); i++) {
    if (i->name() != "h2o" && i->name() != "H2O") {
      M_ += i->molality();
    }
  }
}


void ActivityModel::CalculateActivityCoefficients(
    std::vector<Species>* primary_species,
    std::vector<AqueousEquilibriumComplex>* secondary_species,
    Species* water) {
  const double r1(1.0e0);
  double actw(r1);
  // std::vector<double> gamma;
  // int nsp(primarySpecies->size() + secondarySpecies->size());
  // gamma.resize(nsp, r1);
  // for (auto it = gamma.begin(); it != gamma.end(); ++it) (*it) = r1;

  // Compute activity coefficients
  this->EvaluateVector(*primary_species, *secondary_species, &gamma_, &actw);

  // Set activity coefficients
  int isp(-1);
  for (auto it = primary_species->begin(); it != primary_species->end(); ++it) {
    isp++;
    it->set_act_coef(gamma_[isp]);
    it->update();
  }

  // secondary aqueous complexes
  for (auto it = secondary_species->begin(); it != secondary_species->end(); ++it) {
    isp++;
    it->set_act_coef(gamma_[isp]);
    it->update();
  }

  // Set the water activity
  water->set_act_coef(actw);
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
