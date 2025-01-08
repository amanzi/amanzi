/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  EOS

  EOS for a combination of air and vapor pressure. Mass density is not
  available, not because it can't be calculated, but because it depends
  upon omega. It's not really needed, and if it were, would not fit the
  EOS interface without serious work.
*/

#ifndef AMANZI_EOS_VAPOR_IN_GAS_DENSITY_HH_
#define AMANZI_EOS_VAPOR_IN_GAS_DENSITY_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"

#include "EOS_Density.hh"
#include "Factory.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class VaporInGas_Density : public EOS_Density {
 public:
  VaporInGas_Density(Teuchos::ParameterList& eos_plist);

  virtual double Density(double T, double p) override
  {
    AMANZI_ASSERT(0);
    return 0.0;
  }
  virtual double DDensityDT(double T, double p) override
  {
    AMANZI_ASSERT(0);
    return 0.0;
  }
  virtual double DDensityDp(double T, double p) override
  {
    AMANZI_ASSERT(0);
    return 0.0;
  }

  virtual double MolarDensity(double T, double p) override { return gas_eos_->MolarDensity(T, p); }
  virtual double DMolarDensityDT(double T, double p) override { return gas_eos_->DMolarDensityDT(T, p); }
  virtual double DMolarDensityDp(double T, double p) override { return gas_eos_->DMolarDensityDp(T, p); }

 protected:
  Teuchos::RCP<EOS_Density> gas_eos_;

 private:
  static Utils::RegisteredFactory<EOS_Density, VaporInGas_Density> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
