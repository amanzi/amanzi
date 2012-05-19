/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_SORPTION_ISOTHERM_LINEAR_HH_
#define AMANZI_CHEMISTRY_SORPTION_ISOTHERM_LINEAR_HH_

#include "sorption_isotherm.hh"

// Class for linear isotherm

namespace amanzi {
namespace chemistry {

class SorptionIsothermLinear : public SorptionIsotherm {
 public:
  SorptionIsothermLinear();
  SorptionIsothermLinear(const double KD);
  ~SorptionIsothermLinear();

  void Init(const double KD);
  // returns sorbed concentration
  double Evaluate(const Species& primarySpecies);
  double EvaluateDerivative(const Species& primarySpecies);
  void Display(void) const;

  double KD(void) const { return KD_; }
  void set_KD(const double KD) { KD_ = KD; }

  std::vector<double> GetParameters(void) const;
  void SetParameters(const std::vector<double>& params);

 private:
  // distribution coefficient
  // (currently) units = kg water/m^3 bulk
  double KD_; 

}; // SorptionIsothermLinear

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_SORPTION_ISOTHERM_LINEAR_HH_
