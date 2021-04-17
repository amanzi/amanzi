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

#include "ActivityModel.hh"

namespace Amanzi {
namespace AmanziChemistry {

ActivityModel::ActivityModel()
  : I_(0.0),
    Z_(0.0),
    M_(0.0),
    name_(""),
    num_species_(0) {
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


/* ******************************************************************
* I = 0.5 * sum_i [ m_i * z_i^2 ]
****************************************************************** */
void ActivityModel::CalculateIonicStrength(
    const std::vector<Species>& primary_species,
    const std::vector<AqueousEquilibriumComplex>& secondary_species)
{
  I_ = 0.0;

  // primary species
  for (auto it = primary_species.begin(); it != primary_species.end(); ++it) {
    I_ += it->molality() * it->charge() * it->charge();
  }

  // secondary aqueous complexes
  for (auto it = secondary_species.begin(); it != secondary_species.end(); ++it) {
    I_ += it->molality() * it->charge() * it->charge();
  }

  I_ *= 0.5;
}


/* ******************************************************************
* Z = sum_i [ m_i * |z_i| ]
****************************************************************** */
void ActivityModel::CalculateSumAbsZ(
    const std::vector<Species>& primary_species,
    const std::vector<AqueousEquilibriumComplex>& secondary_species)
{
  Z_ = 0.0;

  // primary species
  for (auto it = primary_species.begin(); it != primary_species.end(); ++it) {
    if (it->name() != "h2o" && it->name() != "H2O") {
      Z_ += it->molality() * std::abs(it->charge());
    }
  }

  // secondary aqueous complexes
  for (auto it = secondary_species.begin(); it != secondary_species.end(); ++it) {
    if (it->name() != "h2o" && it->name() != "H2O") {
      Z_ += it->molality() * std::abs(it->charge());
    }
  }
}


/* ******************************************************************
* TBW
****************************************************************** */
void ActivityModel::CalculateSumC(
    const std::vector<Species>& primary_species,
    const std::vector<AqueousEquilibriumComplex>& secondary_species)
{
  M_ = 0.0;

  // primary species
  for (auto it = primary_species.begin(); it != primary_species.end(); ++it) {
    if (it->name() != "h2o" && it->name() != "H2O") {
      M_ += it->molality();
    }
  }

  // secondary aqueous complexes
  for (auto it = secondary_species.begin(); it != secondary_species.end(); ++it) {
    if (it->name() != "h2o" && it->name() != "H2O") {
      M_ += it->molality();
    }
  }
}


/* ******************************************************************
* TBW
****************************************************************** */
void ActivityModel::CalculateActivityCoefficients(
    std::vector<Species>* primary_species,
    std::vector<AqueousEquilibriumComplex>* secondary_species,
    Species* water)
{
  double actw(1.0);

  // Compute activity coefficients
  this->EvaluateVector(*primary_species, *secondary_species, &gamma_, &actw);

  // Set activity coefficients
  int isp(0);
  for (auto it = primary_species->begin(); it != primary_species->end(); ++it, ++isp) {
    it->set_act_coef(gamma_[isp]);
    it->update();
  }

  // secondary aqueous complexes
  for (auto it = secondary_species->begin(); it != secondary_species->end(); ++it, ++isp) {
    it->set_act_coef(gamma_[isp]);
    it->update();
  }

  // Set the water activity
  water->set_act_coef(actw);
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
