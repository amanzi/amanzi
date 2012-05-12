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
      M_(0.0),
      verbosity_(kSilent),
      name_("") {
}  // end ActivityModel constructor

ActivityModel::~ActivityModel() {
}  // end ActivityModel destructor

void ActivityModel::Setup(
    const ActivityModelParameters& parameters,
    const std::vector<Species>& primary_species,
    const std::vector<AqueousEquilibriumComplex>& secondary_species) {

}  // end Setup()

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
	  if (i->name()!="h2o" && i->name()!="H2O") Z_ += i->molality() * abs(i->charge());
  }

  // secondary aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::const_iterator i = secondarySpecies.begin();
       i != secondarySpecies.end(); i++) {
	  if (i->name()!="h2o" && i->name()!="H2O") Z_ += i->molality() * abs(i->charge());
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
    if (i->name()!="h2o" && i->name()!="H2O") M_ += i->molality();
  }

  // secondary aqueous complexes
  for (std::vector<AqueousEquilibriumComplex>::const_iterator i = secondarySpecies.begin();
       i != secondarySpecies.end(); i++) {
	if (i->name()!="h2o" && i->name()!="H2O") M_ += i->molality();
  }

}  // end CalculateSumAbsZ()

void ActivityModel::CalculateActivityCoefficients(
    std::vector<Species>* primarySpecies,
    std::vector<AqueousEquilibriumComplex>* secondarySpecies,
    Species* water) {
//-----------------------------------------------------
const  double r0(0.0e0), r1(1.0e0);
std::vector<double> gamma;
double actw(r1);
int nsp(primarySpecies->size()+secondarySpecies->size());
gamma.resize(nsp,r1);
for (std::vector<double>::iterator i=gamma.begin(); i!=gamma.end(); i++) (*i)=r1;
//----------------------------------------------------------------------
// Compute activity coefficients
//----------------------------------------------------------------------
this->EvaluateVector (gamma,actw,*primarySpecies,*secondarySpecies);
//----------------------------------------------------------------------
// Set activity coefficients
//----------------------------------------------------------------------
int isp(-1);
for (std::vector<Species>::iterator i = primarySpecies->begin();
       i != primarySpecies->end(); i++) {
	 isp++;
     i->act_coef(gamma[isp]);
     i->update();
}

// secondary aqueous complexes
for (std::vector<AqueousEquilibriumComplex>::iterator i = secondarySpecies->begin();
       i != secondarySpecies->end(); i++) {
      isp++;
	  i->act_coef(gamma[isp]);
	  i->update();
}
// Set the water activity
water->act_coef(actw);

}  // end CalculateActivityCoefficients()

}  // namespace chemistry
}  // namespace amanzi
