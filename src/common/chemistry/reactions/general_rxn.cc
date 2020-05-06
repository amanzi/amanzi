/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for general forward/reverse reaction
*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "general_rxn.hh"
#include "matrix_block.hh"

namespace Amanzi {
namespace AmanziChemistry {

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
}


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
}


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
}


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
}


void GeneralRxn::addContributionToJacobian(
    MatrixBlock* J,
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
        J->AddValue(species_ids_[i], jcomp, stoichiometry_[i]*tempd);
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
        J->AddValue(species_ids_[i], jcomp, stoichiometry_[i]*tempd);
      }
    }  // end columns
  }  // end backward expression
}


void GeneralRxn::display(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  for (unsigned int i = 0; i < species_names_.size(); i++) {
    message << stoichiometry_.at(i) << " " << species_names_.at(i);
    if (i < species_names_.size() - 1) {
      message << " + ";
    }
  }
  message << std::endl;
  message << "        forward_rate = " << std::exp(lnQkf_) << std::endl;
  message << "        backward_rate = " << std::exp(lnQkb_) << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void GeneralRxn::Display(const Teuchos::Ptr<VerboseObject> vo) const {
  // convention for this reaction is that reactants have negative
  // stoichiometries, products have positive stoichiometries....
  // write them in standard chemistry notation by printing -stoich

  std::stringstream message;

  // write the overall reaction
  // reactants:
  message << std::setw(6) << std::fixed << std::setprecision(2);
  for (unsigned int i = 0; i < species_names_.size(); i++) {
    if (stoichiometry_.at(i) < 0) { 
      message << -stoichiometry_.at(i) << " " << species_names_.at(i);
      if (i < forward_species_ids_.size() - 1) {
        message << " + ";
      }
    }
  }
  
  message << " <---> ";
  // products
  for (unsigned int i = 0; i < species_names_.size(); i++) {
    if (stoichiometry_.at(i) > 0) { 
      message << stoichiometry_.at(i) << " " << species_names_.at(i);
      if (i < species_names_.size() - 1) {
        message << " + ";
      }
    }
  }
  message << std::endl;
  message << std::setprecision(6);
  // write the forward rate expression....
  message << std::setw(12) << "    R_f = "
          << std::scientific << this->kf_ << std::fixed;
  if (forward_species_ids_.size() > 0 && this->kf_ > 0.0) {
    message << " * ";
    
    for (unsigned int i = 0; i < forward_species_ids_.size(); i++) {
      message << "a_(" << species_names_[i] << ")^("
              << -stoichiometry_[i] << ")";
      if (i < forward_species_ids_.size() - 1) {
        message << " * ";
      }
    }
  }
  message << std::endl << std::setw(12) << "    R_b = "
          << std::scientific << this->kb_ << std::fixed;
  if (backward_species_ids_.size() > 0 && this->kb_ > 0.0) {
    message << " * ";
    for (unsigned int i = 0; i < backward_species_ids_.size(); i++) {
      message << "a_(" << species_names_[i] << ")^("
              << -stoichiometry_[i] << ")";
      if (i < backward_species_ids_.size() - 1) {
        message << " * ";
      }
    }
  }
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
