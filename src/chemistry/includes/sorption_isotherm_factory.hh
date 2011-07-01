/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_SORP_ISOTHERM_FACTORY_HH_
#define AMANZI_CHEMISTRY_SORP_ISOTHERM_FACTORY_HH_

#include <string>

#include "species.hh"
#include "verbosity.hh"

namespace amanzi {
namespace chemistry {

class SorptionIsotherm;

class SorptionIsothermFactory {
 public:
  SorptionIsothermFactory();
  ~SorptionIsothermFactory();

  SorptionIsotherm* Create(const std::string& model);

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

 protected:

 private:
  Verbosity verbosity_;
};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_SORP_ISOTHERM_FACTORY_HH_
