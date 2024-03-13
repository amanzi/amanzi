/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  Binary diffusion coeffcient, e.g. diffusion of water vapor in air.
*/

#ifndef AMANZI_EOS_VAPOR_IN_GAS_DIFFUSION_HH_
#define AMANZI_EOS_VAPOR_IN_GAS_DIFFUSION_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "Factory.hh"

#include "EOS_Diffusion.hh"

namespace Amanzi {
namespace AmanziEOS {

class VaporInGas_Diffusion : public EOS_Diffusion {
 public:
  explicit VaporInGas_Diffusion(Teuchos::ParameterList& plist);

  virtual double Diffusion(double T, double p);
  virtual double DDiffusionDT(double T, double p);
  virtual double DDiffusionDp(double T, double p);

 private:
  double Tref_; // [K]
  double dref_; // [m^2/s]

  static Utils::RegisteredFactory<EOS_Diffusion, VaporInGas_Diffusion> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
