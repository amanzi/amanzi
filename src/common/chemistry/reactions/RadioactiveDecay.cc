/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry

  Class for radioactive decay of aqueous and sorbed components.
  Does not deal with decay of solid phase.
*/

#include <cmath>
#include <cassert>
#include <sstream>
#include <iostream>
#include <iomanip>

// Amanzi
#include "errors.hh"
#include "StringExt.hh"

// Chemistry
#include "ChemistryUtilities.hh"
#include "MatrixBlock.hh"
#include "RadioactiveDecay.hh"

namespace Amanzi {
namespace AmanziChemistry {

/*************************************************************
* Radioactive decay of aqueous and sorbed components.
****************************************************************** */
RadioactiveDecay::RadioactiveDecay()
  : species_names_(),
    species_ids_(),
    stoichiometry_(),
    rate_constant_(0.0),
    half_life_(0.0),
    rate_(0.0)
{
  ConvertHalfLifeToRateConstant();
}


RadioactiveDecay::RadioactiveDecay(const Teuchos::ParameterList& plist,
                                   const std::map<std::string, int>& name_to_id)
  : rate_constant_(0.0), rate_(0.0)
{
  species_names_.clear();
  species_ids_.clear();
  stoichiometry_.clear();

  half_life_ = plist.get<double>("half life");
  assert(half_life_ > 0.0);

  std::string parent = plist.get<std::string>("reactant");
  int parent_id = name_to_id.at(parent);
  if (parent_id < 0) {
    std::stringstream ss;
    ss << "Unknown parent species '" << parent << "'.\n"
       << "Parent species must be in the primary species list.\n";
    Exceptions::amanzi_throw(Errors::Message(ss.str()));
  }
  species_names_.push_back(parent);
  species_ids_.push_back(parent_id);
  stoichiometry_.push_back(-1.0);

  std::string progeny = plist.get<std::string>("product");
  trim(progeny);

  // NOTE: we allow empty progeny
  if (progeny.size() > 0) {
    double coeff;
    std::string name;
    std::istringstream iss2(progeny);
    while (iss2 >> coeff || !iss2.eof()) {
      iss2 >> name;
      int id2 = name_to_id.at(name);

      if (id2 < 0) {
        std::stringstream ss;
        ss << "Unknown product species '" << progeny << "'.\n"
           << "Product species must be in the primary species list.\n";
        Exceptions::amanzi_throw(Errors::Message(ss.str()));
      }

      species_names_.push_back(name);
      species_ids_.push_back(id2);
      stoichiometry_.push_back(coeff);
    }

    if (species_ids_.size() < 2) {
      std::stringstream ss;
      ss << "Mis-formatted product part of equation '" << progeny << "'.\n";
      Exceptions::amanzi_throw(Errors::Message(ss.str()));
    }
  }

  // we assume that species_names[0] etc is for the parent, any
  // following species are the progeny. The stoichiometry of the
  // parent should be negative!
  assert(species_names_.size() > 0);
  assert(stoichiometry_.at(0) < 0);
  ConvertHalfLifeToRateConstant();
}


/* ******************************************************************
* Solve:
*    C = C_0 * exp(-k*t) for k where C = 0.5*C_0 and t = half life > 0
*    k = -ln(0.5) / half_life
*  The reaction rate constant will be positive (-k is decay) !
****************************************************************** */
void
RadioactiveDecay::ConvertHalfLifeToRateConstant()
{
  rate_constant_ = -std::log(0.5) / half_life_;
}


/* ******************************************************************
*
****************************************************************** */
void
RadioactiveDecay::UpdateRate(const std::vector<double>& total,
                             const std::vector<double>& total_sorbed,
                             const double porosity,
                             const double saturation,
                             const double bulk_volume)
{
  // NOTE: we are working on the totals, not free ion!
  // need total moles of the decaying species:
  //    total * volume_h2o + total_sorbed * volume_bulk
  double volume_h2o = porosity * saturation * bulk_volume * 1000.0; // [L]
  rate_ = volume_h2o * total.at(parent_id());
  if (total_sorbed.size() > 0) { rate_ += bulk_volume * total_sorbed.at(parent_id()); }
  rate_ *= rate_constant();
}


void
RadioactiveDecay::AddContributionToResidual(std::vector<double>* residual)
{
  // Note: rate is < 0, so we add to parent, subtract from
  // progeny. Stoichiometric coefficients should account for this.
  for (int i = 0; i < species_ids_.size(); ++i) {
    int icomp = species_ids_.at(i);
    // this stoichiometry is for the overall reaction
    residual->at(icomp) -= stoichiometry_.at(i) * rate();
  }
}


void
RadioactiveDecay::AddContributionToJacobian(const MatrixBlock& dtotal,
                                            const MatrixBlock& dtotal_sorbed,
                                            const double porosity,
                                            const double saturation,
                                            const double bulk_volume,
                                            MatrixBlock* J)
{
  // NOTE: operating on total [moles/L] not free [moles/kg]
  double volume_h2o = porosity * saturation * bulk_volume * 1000.0; // [L]
  // taking derivative of contribution to residual in row i with respect
  // to species in column j

  // column loop
  int icomp0 = species_ids_.at(0);
  for (int i = 0; i < species_ids_.size(); ++i) {
    int icomp = species_ids_.at(i);
    // row loop
    for (int j = 0; j < J->size(); ++j) {
      double tempd = dtotal(icomp0, j) * volume_h2o;
      if (dtotal_sorbed.size() > 0) { tempd += dtotal_sorbed(icomp0, j) * bulk_volume; }
      tempd *= -rate_constant() * stoichiometry_.at(i);
      J->AddValue(icomp, j, tempd);
    }
  }
}


void
RadioactiveDecay::Display(const Teuchos::Ptr<VerboseObject> vo) const
{
  // convention for this reaction is that reactants have negative
  // stoichiometries, products have positive stoichiometries....
  // write them in standard chemistry notation by printing -stoich

  std::stringstream message;

  // write the overall reaction
  // reactants:
  message << std::setw(6) << std::fixed << std::setprecision(2);
  if (stoichiometry_.at(0) != -1.0) { message << -stoichiometry_.at(0) << " "; }
  message << parent_name() << " --> ";

  // products, note we start at 1 (0=parent)!
  for (unsigned int i = 1; i < species_names_.size(); i++) {
    message << stoichiometry_.at(i) << " " << species_names_.at(i);
    if (i < species_names_.size() - 1) { message << " + "; }
  }
  message << std::endl;
  message << std::setprecision(6);
  // write the rate data
  message << std::setw(20) << "Half Life : ";
  message << std::setw(10) << half_life_ << " [s]\n";
  message << std::setw(20) << " k : " << std::scientific << rate_constant() << std::fixed
          << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

} // namespace AmanziChemistry
} // namespace Amanzi
