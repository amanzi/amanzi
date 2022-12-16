/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry

  Base class for sorption isotherms
*/

#include <string>

#include "SorptionIsotherm.hh"

namespace Amanzi {
namespace AmanziChemistry {

SorptionIsotherm::SorptionIsotherm(const std::string& name, const SorptionIsothermType type)
  : name_(name), isotherm_type_(type)
{}

} // namespace AmanziChemistry
} // namespace Amanzi
