/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "general_rxn.hh"

#include <cmath>

#include <iostream>

#include "block.hh"

namespace amanzi {
namespace chemistry {

GeneralRxn::GeneralRxn() {
  ncomp_ = 0;
  ncomp_forward_ = 0;
  ncomp_backward_ = 0;
  species_names_.clear();
  species_ids_.clear();
  stoichiometry_.clear();
  forward_species_ids_.clear();
  forward_stoichiometry_.clear();
  backward_species_ids_.clear();
  backward_stoichiometry_.clear();
  kf_ = 0.;
  kb_ = 0.;

  lnQkf_ = 0.;
  lnQkb_ = 0.;
}  // end GeneralRxn() constructor

GeneralRxn::GeneralRxn(std::string s) {
  static_cast<void>(s);
}  // end GeneralRxn() constructor

GeneralRxn::GeneralRxn(SpeciesName name,
                       std::vector<SpeciesName>species,
                       std::vector<double>stoichiometries,
                       std::vector<int>species_ids,
                       std::vector<double>forward_stoichiometries,
                       std::vector<int>forward_species_ids,
                       std::vector<double>backward_stoichiometries,
                       std::vector<int>backward_species_ids,
                       double kf, double kb) {
  static_cast<void>(name);
  ncomp_ = species_ids.size();
  ncomp_forward_ = forward_species_ids.size();
  ncomp_backward_ = backward_species_ids.size();

  for (std::vector<SpeciesName>::const_iterator i = species.begin();
       i != species.end(); i++) {
    species_names_.push_back(*i);
  }
  for (std::vector<double>::const_iterator i = stoichiometries.begin();
       i != stoichiometries.end(); i++) {
    stoichiometry_.push_back(*i);
  }
  for (std::vector<int>::const_iterator i = species_ids.begin();
       i != species_ids.end(); i++) {
    species_ids_.push_back(*i);
  }
  // forward species
  for (std::vector<double>::const_iterator i = forward_stoichiometries.begin();
       i != forward_stoichiometries.end(); i++) {
    forward_stoichiometry_.push_back(*i);
  }
  for (std::vector<int>::const_iterator i = forward_species_ids.begin();
       i != forward_species_ids.end(); i++) {
    forward_species_ids_.push_back(*i);
  }
  // backward species
  for (std::vector<double>::const_iterator i = backward_stoichiometries.begin();
       i != backward_stoichiometries.end(); i++) {
    backward_stoichiometry_.push_back(*i);
  }
  for (std::vector<int>::const_iterator i = backward_species_ids.begin();
       i != backward_species_ids.end(); i++) {
    backward_species_ids_.push_back(*i);
  }

  kf_ = kf;
  kb_ = kb;
}  // end GeneralRxn() constructor

GeneralRxn::~GeneralRxn() {
}  // end GeneralRxn() destructor

// temporary location for member functions
void GeneralRxn::update_rates(const std::vector<Species> primarySpecies) {
  // forward rate expression
  lnQkf_ = 0.;
  if (kf_ > 0.) {
    lnQkf_ = std::log(kf_);
    for (int i = 0; i < ncomp_forward_; i++) {
      lnQkf_ += forward_stoichiometry_[i] *
          primarySpecies[ forward_species_ids_[i] ].ln_activity();
    }  // end forward species
  }  // end forward expression

  // backward rate expression
  lnQkb_ = 0.;
  if (kb_ > 0.) {
    lnQkb_ = std::log(kb_);
    for (int i = 0; i < ncomp_backward_; i++) {
      lnQkb_ += backward_stoichiometry_[i] *
          primarySpecies[ backward_species_ids_[i] ].ln_activity();
    }  // end backward species
  }  // end backward expression
}  // end update_rates()

void GeneralRxn::addContributionToResidual(std::vector<double> *residual,
                                           double por_den_sat_vol) {
  // por_den_sat_vol = porosity*water_density*saturation*volume
  double effective_rate = 0.;
  if (kf_ > 0.) {
    effective_rate += std::exp(lnQkf_);
  }
  if (kb_ > 0.) {
    effective_rate -= std::exp(lnQkb_);
  }
  effective_rate *= por_den_sat_vol;

  for (int i = 0; i < ncomp_; i++) {
    int icomp = species_ids_[i];
    // this stoichiometry is for the overall reaction
    (*residual)[icomp] -= stoichiometry_[i] * effective_rate;
  }
}  // end addContributionToResidual()

void GeneralRxn::addContributionToJacobian(
    Block* J,
    const std::vector<Species> primarySpecies,
    double por_den_sat_vol) {

  // taking derivative of contribution to residual in row i with respect
  // to species in column j

  // forward rate expression
  if (kf_ > 0.) {
    // column loop
    for (int j = 0; j < ncomp_forward_; j++) {
      int jcomp = forward_species_ids_[j];
      double tempd = -forward_stoichiometry_[j] *
          std::exp(lnQkf_ - primarySpecies[jcomp].ln_molality()) *
          por_den_sat_vol;
      // row loop
      for (int i = 0; i < ncomp_; i++) {
        J->addValue(species_ids_[i], jcomp, stoichiometry_[i]*tempd);
      }
    }  // end columns
  }  // end forward expression

  // backward rate expression
  if (kb_ > 0.) {
    // column loop
    for (int j = 0; j < ncomp_backward_; j++) {
      int jcomp = backward_species_ids_[j];
      double tempd = backward_stoichiometry_[j] *
          std::exp(lnQkb_ - primarySpecies[jcomp].ln_molality()) *
          por_den_sat_vol;
      // row loop
      for (int i = 0; i < ncomp_; i++) {
        J->addValue(species_ids_[i], jcomp, stoichiometry_[i]*tempd);
      }
    }  // end columns
  }  // end backward expression
}  // end addContributionToJacobian()

void GeneralRxn::display(void) const {
  for (unsigned int i = 0; i < species_names_.size(); i++) {
    std::cout << stoichiometry_[i] << " " << species_names_[i];
    if (i < species_names_.size() - 1) {
      std::cout << " + ";
    }
  }
  std::cout << std::endl;
  std::cout << "        forward_rate = " << std::exp(lnQkf_) << std::endl;
  std::cout << "        backward_rate = " << std::exp(lnQkb_) << std::endl;
}  // end display()

}  // namespace chemistry
}  // namespace amanzi
