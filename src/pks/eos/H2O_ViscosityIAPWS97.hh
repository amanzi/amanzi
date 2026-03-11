/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Viscosity for liquid water consistent with IAPWS97 formulation.
  http://www.iapws.org/relguide/viscosity.html, formulas (10)-(12).
*/

#ifndef AMANZI_EOS_H2O_VISCOSITY_IAPWS97_HH_
#define AMANZI_EOS_H2O_VISCOSITY_IAPWS97_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "Factory.hh"

#include "IAPWS97.hh"
#include "EOS_Viscosity.hh"

namespace Amanzi {
namespace AmanziEOS {

class H2O_ViscosityIAPWS97 : public EOS_Viscosity {
 public:
  explicit H2O_ViscosityIAPWS97(Teuchos::ParameterList& plist);

  virtual double Viscosity(double T, double p);
  virtual double DViscosityDT(double T, double p);
  virtual double DViscosityDp(double T, double p);

 private:
  Teuchos::RCP<IAPWS97> eos_;
  static Utils::RegisteredFactory<EOS_Viscosity, H2O_ViscosityIAPWS97> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
