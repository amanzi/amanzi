/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

#ifndef AMANZI_CHEMISTRY_SORP_ISOTHERM_FACTORY_HH_
#define AMANZI_CHEMISTRY_SORP_ISOTHERM_FACTORY_HH_

#include <string>

#include "species.hh"
#include "string_tokenizer.hh"

namespace Amanzi {
namespace AmanziChemistry {

class SorptionIsotherm;

class SorptionIsothermFactory {
 public:
  SorptionIsothermFactory();
  ~SorptionIsothermFactory() {};

  SorptionIsotherm* Create(const std::string& model, 
                           const StringTokenizer parameters);

  SpeciesId VerifySpeciesName(const SpeciesName species_name,
                              const std::vector<Species>& species) const;

  void set_verbosity(const Verbosity s_verbosity) {
    this->verbosity_ = s_verbosity;
  };
  Verbosity verbosity(void) const {
    return this->verbosity_;
  };

  static const std::string linear;
  static const std::string langmuir;
  static const std::string freundlich;

 private:
  Verbosity verbosity_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
