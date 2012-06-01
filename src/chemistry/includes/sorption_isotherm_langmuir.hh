/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_SORPTION_ISOTHERM_LANGMUIR_HH_
#define AMANZI_CHEMISTRY_SORPTION_ISOTHERM_LANGMUIR_HH_

#include "sorption_isotherm.hh"

// Class for Langmuir isotherm

namespace amanzi {
namespace chemistry {

class SorptionIsothermLangmuir : public SorptionIsotherm {
 public:
  SorptionIsothermLangmuir();
  SorptionIsothermLangmuir(const double K, const double b);
  ~SorptionIsothermLangmuir();

  void Init(const double K, const double b);
  // returns sorbed concentration
  double Evaluate(const Species& primarySpecies);
  double EvaluateDerivative(const Species& primarySpecies);
  void Display(void) const;

  double K(void) const { return K_; }
  void set_K(const double K) { K_ = K; }
  double b(void) const { return b_; }
  void set_b(const double b) { b_ = b; }

  std::vector<double> GetParameters(void) const;
  void SetParameters(const std::vector<double>& params);

private:
  // equilibrium constant or Langmuir adsorption constant
  // units = L water/mol
  double K_; 
  // number of sorption sites (max sorbed concentration)
  // units = mol/m^3 bulk
  double b_; 

}; // SorptionIsothermLangmuir

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_SORPTION_ISOTHERM_LANGMUIR_HH_
