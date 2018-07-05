/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for Langmuir isotherm
*/

#ifndef AMANZI_CHEMISTRY_SORPTION_ISOTHERM_LANGMUIR_HH_
#define AMANZI_CHEMISTRY_SORPTION_ISOTHERM_LANGMUIR_HH_

#include <vector>

#include "sorption_isotherm.hh"

namespace Amanzi {
namespace AmanziChemistry {

class SorptionIsothermLangmuir : public SorptionIsotherm {
 public:
  SorptionIsothermLangmuir();
  SorptionIsothermLangmuir(const double K, const double b);
  ~SorptionIsothermLangmuir() {};

  void Init(const double K, const double b);
  // returns sorbed concentration
  double Evaluate(const Species& primarySpecies);
  double EvaluateDerivative(const Species& primarySpecies);
  void Display(void) const;

  double K(void) const { return K_; }
  void set_K(const double K) { K_ = K; }
  double b(void) const { return b_; }
  void set_b(const double b) { b_ = b; }

  const std::vector<double>& GetParameters(void);
  void SetParameters(const std::vector<double>& params);

private:
  // equilibrium constant or Langmuir adsorption constant
  // units = L water/mol
  double K_; 
  // number of sorption sites (max sorbed concentration)
  // units = mol/m^3 bulk
  double b_;
  std::vector<double> params_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
