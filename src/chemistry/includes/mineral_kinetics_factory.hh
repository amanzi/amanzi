/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_MINERAL_KINETICS_FACTORY_HH_

#define AMANZI_CHEMISTRY_MINERAL_KINETICS_FACTORY_HH_

/*******************************************************************************
 **
 **  File Name: MineralKineticsFactory.h
 **
 **  Description: factory class for building a mineral kinetic rate object
 **
 *******************************************************************************/
#include <vector>
#include <string>

#include "species.hh"
#include "mineral.hh"
#include "verbosity.hh"
#include "string_tokenizer.hh"

namespace amanzi {
namespace chemistry {

class KineticRate;

class MineralKineticsFactory {
 public:
  MineralKineticsFactory(void);
  ~MineralKineticsFactory(void);

  KineticRate* Create(const std::string& rate_type,
                      const StringTokenizer& rate_data,
                      const Mineral& mineral,
                      const SpeciesArray& primary_species);

  SpeciesId VerifyMineralName(const std::string mineral_name,
                              const std::vector<Mineral>& minerals) const;


  void set_verbosity(const Verbosity s_verbosity) {
    this->verbosity_ = s_verbosity;
  };
  Verbosity verbosity(void) const {
    return this->verbosity_;
  };

 protected:

 private:
  Verbosity verbosity_;
  static const std::string kTST;
};

}  // namespace chemistry
}  // namespace amanzi
#endif     /* AMANZI_CHEMISTRY_MINERAL_KINETICS_FACTORY_HH_ */
