/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Constant prescribed viscosity.
*/

#ifndef AMANZI_EOS_VISCOSITY_CONSTANT_HH_
#define AMANZI_EOS_VISCOSITY_CONSTANT_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "Factory.hh"

#include "EOS_Viscosity.hh"

namespace Amanzi {
namespace AmanziEOS {

class ViscosityConstant : public EOS_Viscosity {
 public:
  explicit ViscosityConstant(Teuchos::ParameterList& visc_plist);

  virtual double Viscosity(double T, double p) override { return visc_; }
  virtual double DViscosityDT(double T, double p) override { return 0.0; }
  virtual double DViscosityDp(double T, double p) override { return 0.0; }

 protected:
  virtual void InitializeFromPlist_();

  double visc_;

 private:
  static Utils::RegisteredFactory<EOS_Viscosity, ViscosityConstant> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
