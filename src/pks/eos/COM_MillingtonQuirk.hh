/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  Milliton-Quirk constitutive model for tortuosity.
*/

#ifndef AMANZI_COM_MILLINGTON_QUIRK_HH_
#define AMANZI_COM_MILLINGTON_QUIRK_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "COM_Tortuosity.hh"

namespace Amanzi {
namespace AmanziEOS {

class COM_MillingtonQuirk : public COM_Tortuosity {
 public:
  COM_MillingtonQuirk(Teuchos::ParameterList& plist)
    : COM_Tortuosity(plist), a_(1.0 / 3), b_(7.0 / 3){};
  ~COM_MillingtonQuirk(){};

  virtual double Tortuosity(double phi, double s);
  virtual double DTortuosityDphi(double phi, double s);
  virtual double DTortuosityDs(double phi, double s);

 private:
  double a_, b_;

  static Utils::RegisteredFactory<COM_Tortuosity, COM_MillingtonQuirk> factory_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
