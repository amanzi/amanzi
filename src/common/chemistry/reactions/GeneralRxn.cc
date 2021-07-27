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

#include "GeneralRxn.hh"
#include "MatrixBlock.hh"
#include "ReactionString.hh"

namespace Amanzi {
namespace AmanziChemistry {

GeneralRxn::GeneralRxn() {
  ncomp_ = 0;
  ncomp_forward_ = 0;
  ncomp_backward_ = 0;

  kf_ = 0.0;
  kb_ = 0.0;

  lnQkf_ = 0.0;
  lnQkb_ = 0.0;
}


GeneralRxn::GeneralRxn(const Teuchos::ParameterList& plist,
                       const std::map<std::string, int>& name_to_id)
{
  std::string reactants = plist.get<std::string>("reactants");
  std::string products = plist.get<std::string>("products");

  ParseReaction(reactants, products, &species_names_, &stoichiometry_);

  species_ids_.clear();
  for (auto s = species_names_.begin(); s != species_names_.end(); s++) {
    species_ids_.push_back(name_to_id.at(*s));
  }

  // parse forward rates
  double coeff;
  std::string name;
  forward_species_ids_.clear();
  forward_stoichiometry_.clear();

  std::istringstream iss1(reactants);
  while (iss1 >> coeff || !iss1.eof()) {
    iss1 >> name;
    forward_species_ids_.push_back(name_to_id.at(name));
    forward_stoichiometry_.push_back(coeff);
  }

  // parse backward rates
  backward_species_ids_.clear();
  backward_stoichiometry_.clear();

  std::istringstream iss2(reactants);
  while (iss2 >> coeff || !iss2.eof()) {
    iss2 >> name;
    backward_species_ids_.push_back(name_to_id.at(name));
    backward_stoichiometry_.push_back(coeff);
  }

  kf_ = plist.get<double>("forward rate");
  kb_ = plist.get<double>("backward rate");

  ncomp_ = species_ids_.size();
  ncomp_forward_ = forward_species_ids_.size();
  ncomp_backward_ = backward_species_ids_.size();
}


// temporary location for member functions
void GeneralRxn::UpdateRates(const std::vector<Species> primary_species)
{
  // forward rate expression
  lnQkf_ = 0.0;
  if (kf_ > 0.0) {
    lnQkf_ = std::log(kf_);
    for (int i = 0; i < ncomp_forward_; i++) {
      lnQkf_ += forward_stoichiometry_[i] *
          primary_species[forward_species_ids_[i]].ln_activity();
    }
  }

  // backward rate expression
  lnQkb_ = 0.0;
  if (kb_ > 0.0) {
    lnQkb_ = std::log(kb_);
    for (int i = 0; i < ncomp_backward_; i++) {
      lnQkb_ += backward_stoichiometry_[i] *
          primary_species[backward_species_ids_[i]].ln_activity();
    }
  }
}


void GeneralRxn::AddContributionToResidual(std::vector<double> *residual,
                                           double por_den_sat_vol)
{
  // por_den_sat_vol = porosity*water_density*saturation*volume
  double effective_rate = 0.0;
  if (kf_ > 0.0) {
    effective_rate += std::exp(lnQkf_);
  }
  if (kb_ > 0.0) {
    effective_rate -= std::exp(lnQkb_);
  }
  effective_rate *= por_den_sat_vol;

  for (int i = 0; i < ncomp_; i++) {
    int icomp = species_ids_[i];
    // this stoichiometry is for the overall reaction
    (*residual)[icomp] -= stoichiometry_[i] * effective_rate;
  }
}


void GeneralRxn::AddContributionToJacobian(
    MatrixBlock* J,
    const std::vector<Species> primary_species,
    double por_den_sat_vol)
{
  // taking derivative of contribution to residual in row i with respect
  // to species in column j

  // forward rate expression
  if (kf_ > 0.0) {
    // column loop
    for (int j = 0; j < ncomp_forward_; j++) {
      int jcomp = forward_species_ids_[j];
      double tempd = -forward_stoichiometry_[j] *
          std::exp(lnQkf_ - primary_species[jcomp].ln_molality()) *
          por_den_sat_vol;
      // row loop
      for (int i = 0; i < ncomp_; i++) {
        J->AddValue(species_ids_[i], jcomp, stoichiometry_[i]*tempd);
      }
    }
  }

  // backward rate expression
  if (kb_ > 0.0) {
    // column loop
    for (int j = 0; j < ncomp_backward_; j++) {
      int jcomp = backward_species_ids_[j];
      double tempd = backward_stoichiometry_[j] *
          std::exp(lnQkb_ - primary_species[jcomp].ln_molality()) * por_den_sat_vol;
      // row loop
      for (int i = 0; i < ncomp_; i++) {
        J->AddValue(species_ids_[i], jcomp, stoichiometry_[i]*tempd);
      }
    }
  }
}


void GeneralRxn::Display(const Teuchos::Ptr<VerboseObject> vo) const
{
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
  for (int i = 0; i < species_names_.size(); i++) {
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
    
    for (int i = 0; i < forward_species_ids_.size(); i++) {
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
    for (int i = 0; i < backward_species_ids_.size(); i++) {
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
