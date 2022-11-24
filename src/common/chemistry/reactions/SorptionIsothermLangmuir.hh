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

#include "SorptionIsotherm.hh"

namespace Amanzi {
namespace AmanziChemistry {

class SorptionIsothermLangmuir : public SorptionIsotherm {
 public:
  SorptionIsothermLangmuir();
  SorptionIsothermLangmuir(double K, double b);
  ~SorptionIsothermLangmuir(){};

  void Init(double K, double b);

  // returns sorbed concentration
  virtual double Evaluate(const Species& primary_species) override;
  virtual double EvaluateDerivative(const Species& primary_species) override;

  // setters
  virtual const std::vector<double>& GetParameters() override;
  virtual void SetParameters(const std::vector<double>& params) override;

 private:
  // equilibrium constant or Langmuir adsorption constant [L water/mol]
  double K_;
  // number of sorption sites (max sorbed concentration) [mol/m^3 bulk]
  double b_;
  std::vector<double> params_;
};

} // namespace AmanziChemistry
} // namespace Amanzi
#endif
