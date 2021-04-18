/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Factory class for building a mineral kinetic rate object
*/
 
#include <cstdlib>
#include <string>
#include <sstream>

#include "errors.hh"

#include "KineticRate.hh"
#include "KineticRateFactory.hh"
#include "KineticRateTST.hh"
#include "Mineral.hh"
#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

KineticRate* KineticRateFactory::Create(const Teuchos::ParameterList& plist,
                                        const Mineral& mineral,
                                        const SpeciesArray& primary_species)
{
  std::string model = plist.get<std::string>("rate model");
  double rate = plist.get<double>("rate constant");
  std::string modifiers = plist.get<std::string>("modifiers");

  if (model == "TST") {
    auto kinetic_rate = new KineticRateTST();
    kinetic_rate->Setup(mineral, rate, modifiers, primary_species);
    return kinetic_rate;
  } else {
    std::ostringstream oss;
    oss << "Unknown kinetic rate model: " << model << ", valid names: TST\n";
    Exceptions::amanzi_throw(Errors::Message(oss.str()));
  }

  return NULL;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
