/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for the ideal gas.
*/

#ifndef AMANZI_EOS_DENSITY_IDEAL_GAS_HH_
#define AMANZI_EOS_DENSITY_IDEAL_GAS_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EOS_Density.hh"

namespace Amanzi {
namespace AmanziEOS {
 
// Equation of State model
class EOS_DensityIdealGas : public EOS_Density {
 public:
  explicit EOS_DensityIdealGas(Teuchos::ParameterList& eos_plist);

  virtual double MolarDensity(double T, double p);
  virtual double DMolarDensityDT(double T, double p);
  virtual double DMolarDensityDp(double T, double p);

  virtual double Density(double T, double p) { return MolarDensity(T, p) * M_; }
  virtual double DDensityDT(double T, double p) { return DMolarDensityDT(T, p) * M_; }
  virtual double DDensityDp(double T, double p) { return DMolarDensityDp(T, p) * M_; }

 protected:
  virtual void InitializeFromPlist_();

  double R_;

 private:
  static Utils::RegisteredFactory<EOS_Density, EOS_DensityIdealGas> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
