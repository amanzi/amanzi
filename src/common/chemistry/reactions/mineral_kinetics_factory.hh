/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Factory class for building a mineral kinetic rate object
*/
 
#ifndef AMANZI_CHEMISTRY_MINERAL_KINETICS_FACTORY_HH_
#define AMANZI_CHEMISTRY_MINERAL_KINETICS_FACTORY_HH_

#include <vector>
#include <string>

#include "mineral.hh"
#include "species.hh"

namespace Amanzi {
namespace AmanziChemistry {

class KineticRate;

class MineralKineticsFactory {
 public:
  MineralKineticsFactory() {};
  ~MineralKineticsFactory() {};

  KineticRate* Create(const Teuchos::ParameterList& plist,
                      const Mineral& mineral,
                      const SpeciesArray& primary_species);
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
