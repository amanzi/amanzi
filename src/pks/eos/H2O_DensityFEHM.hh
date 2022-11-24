/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for liquid water for T between 0.001 and 360 C from FEHM manual
*/

#ifndef AMANZI_EOS_H2O_DENSITY_FEHM_HH_
#define AMANZI_EOS_H2O_DENSITY_FEHM_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EOS_Density.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class H2O_DensityFEHM : public EOS_Density {
 public:
  explicit H2O_DensityFEHM(Teuchos::ParameterList& eos_plist);

  virtual double Density(double T, double p) override;
  virtual double DDensityDT(double T, double p) override;
  virtual double DDensityDp(double T, double p) override;

  virtual double MolarDensity(double T, double p) override { return Density(T, p) / M_; }
  virtual double DMolarDensityDT(double T, double p) override { return DDensityDT(T, p) / M_; }
  virtual double DMolarDensityDp(double T, double p) override { return DDensityDp(T, p) / M_; }

 private:
  // constants for water, hard-coded because it would be crazy to try to come
  // up with names for these in a parameter list...
  // -- temperature dependence of density
  double y0_, z0_, y1_, z1_, y2_, z2_, y3_, z3_, y4_, z4_;
  double y5_, z5_, y6_, z6_, y7_, z7_, y8_, z8_, y9_, z9_;
  double T0_;

  static Utils::RegisteredFactory<EOS_Density, H2O_DensityFEHM> factory_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
