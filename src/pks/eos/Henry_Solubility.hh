/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Henry's solubility constant [mol / m^3 / Pa],
     H = kH0 exp[-kH1 (1/T - 1/T0)].
*/

#ifndef AMANZI_EOS_HENRY_SOLUBILITY_HH_
#define AMANZI_EOS_HENRY_SOLUBILITY_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class Henry_Solubility {
 public:
  explicit Henry_Solubility(Teuchos::ParameterList& plist);

  virtual double Solubility(double T);
  virtual double DSolubilityDT(double T);

 private:
  void Init_(Teuchos::ParameterList& plist);

 protected:
  double kH0_, kH1_, T0_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
