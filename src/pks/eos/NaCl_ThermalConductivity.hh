/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  Thermal conductivity for salt.
*/

#ifndef AMANZI_EOS_SALT_THERMAL_CONDUCTIVITY_HH_
#define AMANZI_EOS_SALT_THERMAL_CONDUCTIVITY_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "EOS_ThermalConductivity.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class NaCl_ThermalConductivity : public EOS_ThermalConductivity {
 public:
  explicit NaCl_ThermalConductivity(Teuchos::ParameterList& plist);
  virtual ~NaCl_ThermalConductivity(){};

  virtual double ThermalConductivity(double T, double phi);
  virtual double DThermalConductivityDT(double T, double phi);
  virtual double DThermalConductivityDPhi(double T, double phi);

 protected:
  double kref_, Tref_;
  bool include_phi_, clipping_;

 private:
  static Utils::RegisteredFactory<EOS_ThermalConductivity, NaCl_ThermalConductivity> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
