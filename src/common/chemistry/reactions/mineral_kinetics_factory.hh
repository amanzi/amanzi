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

#include "species.hh"
#include "mineral.hh"
#include "string_tokenizer.hh"

namespace Amanzi {
namespace AmanziChemistry {

class KineticRate;

class MineralKineticsFactory {
 public:
  MineralKineticsFactory(void);
  ~MineralKineticsFactory(void) {};

  KineticRate* Create(const std::string& rate_type,
                      const StringTokenizer& rate_data,
                      const Mineral& mineral,
                      const SpeciesArray& primary_species);

  int VerifyMineralName(const std::string& mineral_name,
                        const std::vector<Mineral>& minerals) const;


  void set_debug(const bool value) {
    this->debug_ = value;
  };
  bool debug(void) const {
    return this->debug_;
  };

 private:
  bool debug_;
  static const std::string kTST;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
