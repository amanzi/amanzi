/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Constant prescribed thermal conductivity.
*/

#ifndef AMANZI_EOS_THERMAL_CONDUCTIVITY_CONSTANT_HH_
#define AMANZI_EOS_THERMAL_CONDUCTIVITY_CONSTANT_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "Factory.hh"

#include "EOS_ThermalConductivity.hh"

namespace Amanzi {
namespace AmanziEOS {

class ThermalConductivityConstant : public EOS_ThermalConductivity {
 public:
  explicit ThermalConductivityConstant(Teuchos::ParameterList& plist)
    : EOS_ThermalConductivity(plist)
  {
    tc_ = plist_.get<double>("reference conductivity");
  }

  virtual double ThermalConductivity(double p, double T, double phi) override { return tc_; }
  virtual double DThermalConductivityDp(double p, double T, double phi) override { return 0.0; }
  virtual double DThermalConductivityDT(double p, double T, double phi) override { return 0.0; }
  virtual double DThermalConductivityDPhi(double p, double T, double phi) override { return 0.0; }

 protected:
  double tc_;

 private:
  static Utils::RegisteredFactory<EOS_ThermalConductivity, ThermalConductivityConstant> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
