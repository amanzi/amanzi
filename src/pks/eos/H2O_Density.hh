/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  EOS for liquid water for T between 0 and 30 C. Need a reference.
*/

#ifndef AMANZI_EOS_H2O_DENSITY_HH_
#define AMANZI_EOS_H2O_DENSITY_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EOS_Density.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of state for water
class H2O_Density : public EOS_Density {
 public:
  H2O_Density(Teuchos::ParameterList& eos_plist);

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
  const double ka_, kb_, kc_, kd_;
  const double kT0_;

  // -- pressure dependence of density
  const double kalpha_, kp0_;

  static Utils::RegisteredFactory<EOS_Density, H2O_Density> factory_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
