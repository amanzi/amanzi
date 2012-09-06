/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "sorption_isotherm.hh"

#include <iostream>
#include <iomanip>

namespace amanzi {
namespace chemistry {

SorptionIsotherm::SorptionIsotherm(const std::string name,
                                   const SorptionIsothermType type)
    : name_(name),
      isotherm_type_(type) {

}  // end SorptionIsotherm() constructor

SorptionIsotherm::~SorptionIsotherm() {
}  // end SorptionIsotherm() destructor

}  // namespace chemistry
}  // namespace amanzi
