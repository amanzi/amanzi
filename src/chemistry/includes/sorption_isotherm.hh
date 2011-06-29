/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_SORPTION_ISOTHERM_HH_
#define AMANZI_CHEMISTRY_SORPTION_ISOTHERM_HH_

// Base class for sorption isotherms

#include<species.hh>

namespace amanzi {
namespace chemistry {

class SorptionIsotherm {
 public:
  SorptionIsotherm();
  virtual~SorptionIsotherm();

  virtual double Evaluate(const Species& primarySpecies) = 0;
  virtual double EvaluateDerivative(const Species& primarySpecies) = 0;

}; // SorptionIsotherm

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_SORPTION_ISOTHERM_HH_
