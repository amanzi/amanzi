/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  EOS for the ideal gas.
*/

#ifndef AMANZI_EOS_IDEAL_GAS_DENSITY_HH_
#define AMANZI_EOS_IDEAL_GAS_DENSITY_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EOS_Density.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class IdealGas_Density : public EOS_Density {
 public:
  explicit IdealGas_Density(Teuchos::ParameterList& eos_plist);

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
  static Utils::RegisteredFactory<EOS_Density, IdealGas_Density> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
