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

#include "dbc.hh"

#include "GeneralRxn.hh"
#include "MatrixBlock.hh"
#include "ReactionString.hh"

namespace Amanzi {
namespace AmanziChemistry {

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
  species_ids_f_.clear();
  stoichiometry_f_.clear();
  orders_oc_.clear();

  std::istringstream iss1(reactants);
  while (iss1 >> coeff || !iss1.eof()) {
    iss1 >> name;
    species_ids_f_.push_back(name_to_id.at(name));
    stoichiometry_f_.push_back(coeff);
    orders_oc_.push_back(coeff);
  }

  // parse backward rates
  species_ids_b_.clear();
  stoichiometry_b_.clear();

  std::istringstream iss2(products);
  while (iss2 >> coeff || !iss2.eof()) {
    iss2 >> name;
    species_ids_b_.push_back(name_to_id.at(name));
    stoichiometry_b_.push_back(coeff);
    orders_oc_.push_back(coeff);
  }

  kf_ = plist.get<double>("forward rate");
  kb_ = plist.get<double>("backward rate");

  lnkf_ = std::log(kf_);
  lnkb_ = (kb_ > 0.0) ? std::log(kb_) : 0.0;

  ncomp_ = species_ids_.size();
  ncomp_forward_ = species_ids_f_.size();
  ncomp_backward_ = species_ids_b_.size();

  // fractional order in C
  flag_oc_ = false; 
  if (plist.isParameter("reaction orders (reactants/products)")) {
    orders_oc_ = plist.get<Teuchos::Array<double>>("reaction orders (reactants/products)").toVector();
    flag_oc_ = true; 
    AMANZI_ASSERT(orders_oc_.size() == ncomp_);
  }
}


// temporary location for member functions
void GeneralRxn::UpdateRates(const std::vector<Species> primary_species)
{
  // forward rate expression
  lnQkf_ = 0.0;
  lnOcf_ = 0.0;
  if (kf_ > 0.0) {
    for (int i = 0; i < ncomp_forward_; i++) {
      lnQkf_ += stoichiometry_f_[i] * primary_species[species_ids_f_[i]].ln_activity();
    }
    lnOcf_ += lnQkf_;
    lnQkf_ += lnkf_;
  }

  // backward rate expression
  lnQkb_ = 0.0;
  lnOcb_ = 0.0;
  if (kb_ > 0.0) {
    for (int i = 0; i < ncomp_backward_; i++) {
      lnQkb_ += stoichiometry_b_[i] * primary_species[species_ids_b_[i]].ln_activity();
    }
    lnOcb_ += lnQkb_;
    lnQkb_ += lnkb_;
  }

  // fractinal order in C (primary concentrations)
  if (flag_oc_) {
    lnOcf_ = 0.0;
    for (int i = 0; i < ncomp_forward_; i++) {
      lnOcf_ += orders_oc_[i] * primary_species[species_ids_f_[i]].ln_activity();
    }

    lnOcb_ = 0.0;
    for (int i = 0; i < ncomp_backward_; i++) {
      int k = ncomp_forward_ + i;
      lnOcb_ += orders_oc_[k] * primary_species[species_ids_b_[i]].ln_activity();
    }
  }
}


void GeneralRxn::AddContributionToResidual(std::vector<double> *residual,
                                           double por_den_sat_vol)
{
  // por_den_sat_vol = porosity*water_density*saturation*volume
  double effective_rate = 0.0;
  if (kf_ > 0.0) {
    effective_rate += kf_ * std::exp(lnOcf_ + lnOcb_);
  }
  if (kb_ > 0.0) {
    effective_rate -= kf_ * std::exp(lnOcf_ + lnOcb_ + lnQkb_ - lnQkf_);
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
      int jcomp = species_ids_f_[j];
      // double tempd = -stoichiometry_[j] *
      double tempd = orders_oc_[j] *
          std::exp(lnQkf_ - primary_species[jcomp].ln_molality()) * por_den_sat_vol;
      // row loop
      for (int i = 0; i < ncomp_; i++) {
        J->AddValue(species_ids_[i], jcomp, stoichiometry_[i] * tempd);
      }
    }
  }

  // backward rate expression
  if (kb_ > 0.0) {
    // column loop
    for (int j = 0; j < ncomp_backward_; j++) {
      int jcomp = species_ids_b_[j];
      double tempd = orders_oc_[ncomp_forward_ + j] *
          std::exp(lnQkb_ - primary_species[jcomp].ln_molality()) * por_den_sat_vol;
      // row loop
      for (int i = 0; i < ncomp_; i++) {
        J->AddValue(species_ids_[i], jcomp, stoichiometry_[i] * tempd);
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
      if (i < species_ids_f_.size() - 1) {
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
  if (species_ids_f_.size() > 0 && this->kf_ > 0.0) {
    message << " * ";
    
    for (int i = 0; i < species_ids_f_.size(); i++) {
      message << "a_(" << species_names_[i] << ")^(" << orders_oc_[i] << ")";
      if (i < species_ids_f_.size() - 1) {
        message << " * ";
      }
    }
  }
  message << std::endl << std::setw(12) << "    R_b = "
          << std::scientific << this->kb_ << std::fixed;
  if (species_ids_b_.size() > 0 && this->kb_ > 0.0) {
    message << " * ";
    for (int i = 0; i < species_ids_b_.size(); i++) {
      message << "a_(" << species_names_[i] << ")^(" << orders_oc_[ncomp_forward_ + i] << ")";
      if (i < species_ids_b_.size() - 1) {
        message << " * ";
      }
    }
  }
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
