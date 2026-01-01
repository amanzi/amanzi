/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  EOS for liquid water based on the International Association for the 
  Properties of Water and Steam (IAPWS), Industrial Formulation 1997.
*/

#ifndef AMANZI_EOS_H2O_DENSITY_IAPWS97_HH_
#define AMANZI_EOS_H2O_DENSITY_IAPWS97_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "IAPWS97.hh"
#include "EOS_Density.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class H2O_DensityIAPWS97 : public EOS_Density {
 public:
  explicit H2O_DensityIAPWS97(Teuchos::ParameterList& eos_plist);

  virtual double Density(double T, double p) override;
  virtual double DDensityDT(double T, double p) override;
  virtual double DDensityDp(double T, double p) override;

  virtual double MolarDensity(double T, double p) override { return Density(T, p) / M_; }
  virtual double DMolarDensityDT(double T, double p) override { return DDensityDT(T, p) / M_; }
  virtual double DMolarDensityDp(double T, double p) override { return DDensityDp(T, p) / M_; }

 private:
  Teuchos::RCP<IAPWS97> eos_;
  static Utils::RegisteredFactory<EOS_Density, H2O_DensityIAPWS97> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
