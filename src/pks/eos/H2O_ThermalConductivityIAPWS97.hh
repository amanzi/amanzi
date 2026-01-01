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
  http://www.iapws.org/relguide/ThCond.html
*/

#ifndef AMANZI_EOS_H2O_THERMAL_CONDUCTIVITY_IAPWS97_HH_
#define AMANZI_EOS_H2O_THERMAL_CONDUCTIVITY_IAPWS97_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "H2O_ThermalConductivity.hh"
#include "IAPWS97.hh"

namespace Amanzi {
namespace AmanziEOS {

class H2O_ThermalConductivityIAPWS97 : public EOS_ThermalConductivity {
 public:
  explicit H2O_ThermalConductivityIAPWS97(Teuchos::ParameterList& plist);
  virtual ~H2O_ThermalConductivityIAPWS97() {};

  virtual double ThermalConductivity(double p, double T, double phi);
  virtual double DThermalConductivityDp(double p, double T, double phi);
  virtual double DThermalConductivityDT(double p, double T, double phi);
  virtual double DThermalConductivityDPhi(double p, double T, double phi)
  {
    AMANZI_ASSERT(false);
    return 0.0;
  }

 private:
  Teuchos::RCP<IAPWS97> eos_;
  static Utils::RegisteredFactory<EOS_ThermalConductivity, H2O_ThermalConductivityIAPWS97> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
