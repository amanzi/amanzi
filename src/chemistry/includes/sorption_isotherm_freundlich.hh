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

  void Init(const double KD, const double n);
  // returns sorbed concentration
  double Evaluate(const Species& primarySpecies);
  double EvaluateDerivative(const Species& primarySpecies);
  void Display(void) const;

  double KD(void) const { return KD_; }
  void set_KD(const double KD) { KD_ = KD; }
  double n(void) const { return 1./one_over_n_; }
  void set_n(const double n) { one_over_n_ = 1./n; }
  double one_over_n(void) const { return one_over_n_; }
  void set_one_over_n(const double one_over_n) { one_over_n_ = one_over_n; }

 private:
  double KD_; // distribution coefficient
  double one_over_n_; // chemical-specific constant

}; // SorptionIsothermFreundlich

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_SORPTION_ISOTHERM_FREUNDLICH_HH_
