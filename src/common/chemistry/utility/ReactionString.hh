/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Supporting non-member functions
*/

#ifndef AMANZI_CHEMISTRY_BEAKER_UTILITY_HH_
#define AMANZI_CHEMISTRY_BEAKER_UTILITY_HH_

#include <string>
#include <vector>

namespace Amanzi {
namespace AmanziChemistry {

void ParseReaction(const std::string& reactants,
                   const std::string& products,
                   std::vector<std::string>* species,
                   std::vector<double>* stoichiometries);

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
