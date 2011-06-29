/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_SORP_ISOTHERM_FACTORY_HH_
#define AMANZI_CHEMISTRY_SORP_ISOTHERM_FACTORY_HH_

#include <string>

namespace amanzi {
namespace chemistry {

class SorptionIsotherm;

class SorptionIsothermFactory {
 public:
  SorptionIsothermFactory();
  ~SorptionIsothermFactory();

  SorptionIsotherm* Create(const std::string& model);

  static const std::string linear;
  static const std::string langmuir;
  static const std::string freundlich;

 protected:

 private:
};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_SORP_ISOTHERM_FACTORY_HH_
