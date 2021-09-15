/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

#ifndef AMANZI_CHEMISTRY_SORP_ISOTHERM_FACTORY_HH_
#define AMANZI_CHEMISTRY_SORP_ISOTHERM_FACTORY_HH_

#include <memory>
#include <string>

#include "Teuchos_ParameterList.hpp"

#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

class SorptionIsotherm;

class SorptionIsothermFactory {
 public:
  SorptionIsothermFactory() {};
  ~SorptionIsothermFactory() {};

  std::shared_ptr<SorptionIsotherm> Create(const Teuchos::ParameterList& plist);

  int VerifySpeciesName(const std::string& species_name,
                        const std::vector<Species>& species) const;

  static const std::string linear;
  static const std::string langmuir;
  static const std::string freundlich;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
