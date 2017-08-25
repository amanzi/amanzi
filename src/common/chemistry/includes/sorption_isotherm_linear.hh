/* -*-  mode: c++; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_SORPTION_ISOTHERM_LINEAR_HH_
#define AMANZI_CHEMISTRY_SORPTION_ISOTHERM_LINEAR_HH_

#include "sorption_isotherm.hh"

// Class for linear isotherm

namespace Amanzi {
namespace AmanziChemistry {

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

  const std::vector<double>& GetParameters(void);
  void SetParameters(const std::vector<double>& params);

 private:
  // distribution coefficient
  // (currently) units = kg water/m^3 bulk
  double KD_; 
  std::vector<double> params_;

}; // SorptionIsothermLinear

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif  // AMANZI_CHEMISTRY_SORPTION_ISOTHERM_LINEAR_HH_
