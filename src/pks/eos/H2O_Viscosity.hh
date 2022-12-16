/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Viscosity for liquid water.
*/

#ifndef AMANZI_EOS_H2O_VISCOSITY_HH_
#define AMANZI_EOS_H2O_VISCOSITY_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "Factory.hh"

#include "EOS_Viscosity.hh"

namespace Amanzi {
namespace AmanziEOS {

class H2O_Viscosity : public EOS_Viscosity {
 public:
  explicit H2O_Viscosity(Teuchos::ParameterList& eos_plist);

  virtual double Viscosity(double T, double p);
  virtual double DViscosityDT(double T, double p);
  virtual double DViscosityDp(double T, double p);

 protected:
  // constants for water, hard-coded
  const double kav1_, kbv1_, kcv1_;
  const double kbv2_, kcv2_, kT1_;

 private:
  static Utils::RegisteredFactory<EOS_Viscosity, H2O_Viscosity> factory_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
