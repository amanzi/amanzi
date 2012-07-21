/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_SORPTION_ISOTHERM_FREUNDLICH_HH_
#define AMANZI_CHEMISTRY_SORPTION_ISOTHERM_FREUNDLICH_HH_

#include "sorption_isotherm.hh"

// Class for Freundlich isotherm

namespace amanzi {
namespace chemistry {

class SorptionIsothermFreundlich : public SorptionIsotherm {
 public:
  SorptionIsothermFreundlich();
  SorptionIsothermFreundlich(const double KD, const double n);
  ~SorptionIsothermFreundlich();

  // returns sorbed concentration
  double Evaluate(const Species& primarySpecies);
  double EvaluateDerivative(const Species& primarySpecies);
  void Display(void) const;

  double KD(void) const { return KD_; }
  void set_KD(const double KD) { KD_ = KD; }
  double n(void) const { return n_; }
  void set_n(const double n) { n_ = n; }

  std::vector<double> GetParameters(void) const;
  void SetParameters(const std::vector<double>& params);

 private:
  double KD_; // distribution coefficient
  double n_; // chemical-specific constant

}; // SorptionIsothermFreundlich

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_SORPTION_ISOTHERM_FREUNDLICH_HH_
