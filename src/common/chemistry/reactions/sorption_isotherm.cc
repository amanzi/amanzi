/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Base class for sorption isotherms
*/

#include <iostream>
#include <iomanip>

#include "sorption_isotherm.hh"

namespace Amanzi {
namespace AmanziChemistry {

SorptionIsotherm::SorptionIsotherm(const std::string& name,
                                   const SorptionIsothermType type)
    : name_(name),
      isotherm_type_(type) {
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
