/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "sorption_isotherm_rxn.hh"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "chemistry_exception.hh"
#include "block.hh"

#include "exceptions.hh"

namespace amanzi {
namespace chemistry {

SorptionIsothermRxn::SorptionIsothermRxn() {
}

SorptionIsothermRxn::~SorptionIsothermRxn() {
}

void SorptionIsothermRxn::Update(const std::vector<Species>& primarySpecies) {
  sorbed_concentration_= (*isotherm_).Evaluate(primarySpecies.at(species_id_));
}  // end Update()

void SorptionIsothermRxn::AddContributionToTotal(std::vector<double> *total) {
  (*total)[species_id_] += sorbed_concentration_;
}  // end AddContributionToTotal()

void SorptionIsothermRxn::AddContributionToDTotal(
    const std::vector<Species>& primarySpecies,
    Block* dtotal) {

}  // end AddContributionToDTotal()

void SorptionIsothermRxn::Display(void) const {
}  // end Display()

}  // namespace chemistry
}  // namespace amanzi
