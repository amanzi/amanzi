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


void ActivityModel::ResizeGamma(const int size) {
  set_num_species(size);
  gamma_.resize(num_species(), 1.0);  
}

void ActivityModel::CalculateIonicStrength(
    const std::vector<Species>& primarySpecies,
    const std::vector<AqueousEquilibriumComplex>& secondarySpecies) {
  // I = 0.5 * sum_i(m_i*z_i^2)
  I_ = 0.0;

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
    std::vector<Species>* primarySpecies,
    std::vector<AqueousEquilibriumComplex>* secondarySpecies,
    Species* water) {
  const double r1(1.0e0);
  double actw(r1);
  //std::vector<double> gamma;
  // int nsp(primarySpecies->size() + secondarySpecies->size());
  // gamma.resize(nsp, r1);
  // for (std::vector<double>::iterator i = gamma.begin(); i != gamma.end(); i++) {
  //   (*i) = r1;
  // }

  // Compute activity coefficients
  this->EvaluateVector(*primarySpecies, *secondarySpecies, &gamma_, &actw);

  // Set activity coefficients
  int isp(-1);
  for (std::vector<Species>::iterator i = primarySpecies->begin();
       i != primarySpecies->end(); ++i) {
    isp++;
    i->act_coef(gamma_[isp]);
    i->update();
  }

  // secondary aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::iterator i = secondarySpecies->begin();
       i != secondarySpecies->end(); ++i) {
    isp++;
    i->act_coef(gamma_[isp]);
    i->update();
  }
  // Set the water activity
  water->act_coef(actw);
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
