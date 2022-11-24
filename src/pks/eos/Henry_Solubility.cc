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

#include "Henry_Solubility.hh"

namespace Amanzi {
namespace AmanziEOS {

Henry_Solubility::Henry_Solubility(Teuchos::ParameterList& plist)
  : kH0_(7.698e-6), // hydrogen
    kH1_(0.0),
    T0_(298.15)
{
  Init_(plist);
};


double
Henry_Solubility::Solubility(double T)
{
  return (kH1_ == 0.0) ? kH0_ : kH0_ * std::exp(-kH1_ * (1.0 / T - 1.0 / T0_));
};


double
Henry_Solubility::DSolubilityDT(double T)
{
  return (kH1_ == 0.0) ? kH0_ : Solubility(T) * kH1_ / (T * T);
};


/* *******************************************************************
* Solute specific initilazation
******************************************************************* */
void
Henry_Solubility::Init_(Teuchos::ParameterList& plist)
{
  kH0_ = plist.get<double>("reference solubility");
  kH1_ = plist.get<double>("van't Hoff factor", 0.0);
}

} // namespace AmanziEOS
} // namespace Amanzi
