/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "activity_model.hh"

#include <cmath>
#include <cstdlib>
#include <math.h>

#include <vector>
#include <iostream>

namespace amanzi {
namespace chemistry {

ActivityModel::ActivityModel()
    : I_(0.0),
      Z_(0.0),
      M_(0.0) {
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

void ActivityModel::CalculateSumAbsZ(
    const std::vector<Species>& primarySpecies,
    const std::vector<AqueousEquilibriumComplex>& secondarySpecies) {
  // Z = sum_i(m_i*abs(z_i))
  Z_ = 0.0e0;

  // primary species
  for (std::vector<Species>::const_iterator i = primarySpecies.begin();
       i != primarySpecies.end(); i++) {
    Z_ += i->molality() * abs(i->charge());
  }

  // secondary aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::const_iterator i = secondarySpecies.begin();
       i != secondarySpecies.end(); i++) {
    Z_ += i->molality() * abs(i->charge());
  }

}  // end CalculateSumAbsZ()

void ActivityModel::CalculateSumC(
    const std::vector<Species>& primarySpecies,
    const std::vector<AqueousEquilibriumComplex>& secondarySpecies) {
  // Z = sum_i(m_i*abs(z_i))
  M_ = 0.0e0;

  // primary species
  for (std::vector<Species>::const_iterator i = primarySpecies.begin();
       i != primarySpecies.end(); i++) {
    if (i->name()!="h2o")
	   M_ += i->molality();
  }

  // secondary aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::const_iterator i = secondarySpecies.begin();
       i != secondarySpecies.end(); i++) {
	  if (i->name()!="h2o")
	  	   M_ += i->molality();
  }

}  // end CalculateSumAbsZ()

void ActivityModel::CalculateActivityCoefficients(
    std::vector<Species>* primarySpecies,
    std::vector<AqueousEquilibriumComplex>* secondarySpecies) {
//-----------------------------------------------------
const  double r0(0.0e0), r1(1.0e0);
std::vector<double> gamma;
int nsp(primarySpecies->size()+secondarySpecies->size());
gamma.resize(nsp,r1);
for (std::vector<double>::iterator i=gamma.begin(); i!=gamma.end(); i++) (*i)=r1;
//----------------------------------------------------------------------
// Compute activity coefficients
//----------------------------------------------------------------------
this->EvaluateVector (gamma,*primarySpecies,*secondarySpecies);
//----------------------------------------------------------------------
// Set activity coefficients
//----------------------------------------------------------------------
int isp(-1);
for (std::vector<Species>::iterator i = primarySpecies->begin();
       i != primarySpecies->end(); i++) {
	 isp++;
     i->act_coef(gamma[isp]);
}

// secondary aqueous complexes
for (std::vector<AqueousEquilibriumComplex>::iterator i = secondarySpecies->begin();
       i != secondarySpecies->end(); i++) {
      isp++;
	  i->act_coef(gamma[isp]);
}
}  // end CalculateActivityCoefficients()

}  // namespace chemistry
}  // namespace amanzi
