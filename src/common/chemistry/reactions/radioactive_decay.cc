/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for radioactive decay of aqueous and sorbed components.
  Does not deal with decay of solid phase.
*/
 
#include <cmath>
#include <cassert>
#include <sstream>
#include <iostream>
#include <iomanip>

// TPLs
#include "boost/algorithm/string.hpp"

// Chemistry
#include "chemistry_utilities.hh"
#include "chemistry_exception.hh"
#include "matrix_block.hh"
#include "radioactive_decay.hh"

namespace Amanzi {
namespace AmanziChemistry {

/*
**  Simple class for radioactive decay of aqueous and sorbed
**  components. This class is not general enough to handle the decay
**  of elements that are incorporated into minerals.
*/
RadioactiveDecay::RadioactiveDecay()
    : species_names_(),
      species_ids_(),
      stoichiometry_(),
      rate_constant_(0.0),
      half_life_user_(1.0),
      half_life_units_("seconds"),
      half_life_seconds_(0.0),
      rate_(0.0) {
  ConvertHalfLifeUnits();
  ConvertHalfLifeToRateConstant();
}  // end RadioactiveDecay() constructor


RadioactiveDecay::RadioactiveDecay(const std::vector<SpeciesName> species_names,
                                   const std::vector<int> species_ids,
                                   const std::vector<double> stoichiometries,
                                   const double half_life,
                                   const std::string half_life_units)
    : species_names_(species_names),
      species_ids_(species_ids),
      stoichiometry_(stoichiometries),      
      rate_constant_(0.0),
      half_life_user_(half_life),
      half_life_units_(half_life_units),
      half_life_seconds_(0.0),
      rate_(0.0) {
  // we assume that species_names[0] etc is for the parent, any
  // following species are the progeny. The stoichiometry of the
  // parent should be negative!
  assert(species_names_.size() > 0);
  assert(species_ids_.size() > 0);
  assert(stoichiometry_.size() > 0);
  assert(stoichiometry_.at(0) < 0);
  ConvertHalfLifeUnits();
  ConvertHalfLifeToRateConstant();
}


void RadioactiveDecay::ConvertHalfLifeUnits(void) {
  double conversion = 1.0;
  std::string units = half_life_units_;
  utilities::RemoveLeadingAndTrailingWhitespace(&units);
  if (boost::iequals(units, "years") || boost::iequals(units, "y")) {
    conversion = 365.0 * 24.0 * 60.0 * 60.0; 
  } else if (boost::iequals(units, "days") || boost::iequals(units, "d")) {
    conversion = 24.0 * 60.0 * 60.0;
  } else if (boost::iequals(units, "hours") || boost::iequals(units, "h")) {
    conversion = 60.0 * 60.0;
  } else if (boost::iequals(units, "minutes") || boost::iequals(units, "m")) {
    conversion = 60.0;
  } else if (boost::iequals(units, "seconds") || boost::iequals(units, "s")) {
    conversion = 1.0;
  } else {
    std::stringstream message;
    message << "ERROR: RadioactiveDecay::ConvertHalfLifeUnits(" << parent_name() << "): \n"
            << "Unknown half life units '" << half_life_units_ << "'.\n"
            << "Valid units are: 'years', 'days', 'hours', 'minutes', 'seconds'.\n";
    Exceptions::amanzi_throw(ChemistryInvalidInput(message.str()));
  }

  half_life_seconds_ = half_life_user_ * conversion;
}


void RadioactiveDecay::ConvertHalfLifeToRateConstant(void) {
  /*
  **  Solve:
  **    C = C_0 * exp(-k*t) for k where C = 0.5*C_0 and t = half life > 0
  **    k = -ln(0.5) / half_life
  **  The reaction rate constant will be positive (-k is decay)!
  */
  rate_constant_ = -std::log(0.5) / half_life_seconds_;
  assert(rate_constant_ > 0.0);
}


// temporary location for member functions
void RadioactiveDecay::UpdateRate(const std::vector<double>& total,
                                  const std::vector<double>& total_sorbed,
                                  const double porosity,
                                  const double saturation,
                                  const double bulk_volume) {
  // NOTE: we are working on the totals, not free ion!
  // need total moles of the decaying species:
  //    total * volume_h2o + total_sorbed * volume_bulk
  double volume_h2o = porosity * saturation * bulk_volume * 1000.0; // [L]
  rate_ = volume_h2o * total.at(parent_id());
  if (total_sorbed.size() > 0) {
    rate_ += bulk_volume * total_sorbed.at(parent_id());
  }
  rate_ *= rate_constant();
}
 

void RadioactiveDecay::AddContributionToResidual(std::vector<double> *residual) {
  // Note: rate is < 0, so we add to parent, subtract from
  // progeny. Stoichiometric coefficients should account for this.
  for (int i = 0; i < species_ids_.size(); ++i) {
    int icomp = species_ids_.at(i);
    // this stoichiometry is for the overall reaction
    residual->at(icomp) -= stoichiometry_.at(i) * rate();
  }
}


void RadioactiveDecay::AddContributionToJacobian(
    const MatrixBlock& dtotal, const MatrixBlock& dtotal_sorbed,
    const double porosity, const double saturation, const double bulk_volume,
    MatrixBlock* J) {

  // NOTE: operating on total [moles/L] not free [moles/kg]
  double volume_h2o = porosity * saturation * bulk_volume * 1000.0; // [L]
  // taking derivative of contribution to residual in row i with respect
  // to species in column j

  // column loop
  for (int i = 0; i < species_ids_.size(); ++i) {
    int icomp = species_ids_.at(i);
    // row loop
    for (int j = 0; j < J->size(); ++j) {
      double tempd = dtotal(icomp, j) * volume_h2o;
      if (dtotal_sorbed.size() > 0) {
        tempd += dtotal_sorbed(icomp, j) * bulk_volume;
      }
      tempd *= -rate_constant() * stoichiometry_.at(i);
      J->AddValue(icomp, j, tempd);
    }
  }  // end columns
}


void RadioactiveDecay::Display(const Teuchos::Ptr<VerboseObject> vo) const {
  // convention for this reaction is that reactants have negative
  // stoichiometries, products have positive stoichiometries....
  // write them in standard chemistry notation by printing -stoich

  std::stringstream message;

  // write the overall reaction
  // reactants:
  message << std::setw(6) << std::fixed << std::setprecision(2);
  if (stoichiometry_.at(0) != -1.0) {
    message << -stoichiometry_.at(0) << " ";
  }
  message << parent_name() << " --> ";
  // products, note we start at 1 (0=parent)!
  for (unsigned int i = 1; i < species_names_.size(); i++) {
    message << stoichiometry_.at(i) << " " << species_names_.at(i);
    if (i < species_names_.size() - 1) {
      message << " + ";
    }
  }
  message << std::endl;
  message << std::setprecision(6);
  // write the rate data
  message << std::setw(20) << "Half Life : ";
  message << std::setw(10) << half_life_user_ << " [" << half_life_units_ << "]    ";
  message << std::setw(10) << std::scientific << half_life_seconds_ << std::fixed << " [seconds]" << std::endl;
  message << std::setw(20) << " k : " << std::scientific << rate_constant() << std::fixed << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());

}

}  // namespace AmanziChemistry
}  // namespace Amanzi
