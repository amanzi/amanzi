/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_SORPTION_ISOTHERM_HH_
#define AMANZI_CHEMISTRY_SORPTION_ISOTHERM_HH_

// Base class for sorption isotherms

#include <string>

#include<species.hh>

namespace amanzi {
namespace chemistry {

class SorptionIsotherm {
 public:
  SorptionIsotherm(const std::string name);
  virtual~SorptionIsotherm();

  virtual double Evaluate(const Species& primarySpecies) = 0;
  virtual double EvaluateDerivative(const Species& primarySpecies) = 0;

  virtual void Display(void) const = 0;

  virtual std::vector<double> GetParameters(void) const = 0;
  virtual void SetParameters(const std::vector<double>& params) = 0;

  std::string name(void) const {
    return name_;
  }

 protected:

 private:
  std::string name_;

}; // SorptionIsotherm

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_SORPTION_ISOTHERM_HH_
