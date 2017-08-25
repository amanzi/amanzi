/* -*-  mode: c++; indent-tabs-mode: nil -*- */

#include "sorption_isotherm.hh"

#include <iostream>
#include <iomanip>

namespace Amanzi {
namespace AmanziChemistry {

SorptionIsotherm::SorptionIsotherm(const std::string name,
                                   const SorptionIsothermType type)
    : name_(name),
      isotherm_type_(type) {

}  // end SorptionIsotherm() constructor

SorptionIsotherm::~SorptionIsotherm() {
}  // end SorptionIsotherm() destructor

}  // namespace AmanziChemistry
}  // namespace Amanzi
