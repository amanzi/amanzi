/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  Constitutive model for tortuosity \tau(\phi, s) as the function
  of porosity and liquid saturation.
*/

#ifndef AMANZI_COM_TORTUOSITY_HH_
#define AMANZI_COM_TORTUOSITY_HH_

namespace Amanzi {
namespace AmanziEOS {

#include "Teuchos_ParameterList.hpp"

class COM_Tortuosity {
 public:
  COM_Tortuosity(Teuchos::ParameterList& plist) : plist_(plist){};
  virtual ~COM_Tortuosity(){};

  virtual double Tortuosity(double phi, double s) = 0;
  virtual double DTortuosityDphi(double phi, double s) = 0;
  virtual double DTortuosityDs(double phi, double s) = 0;

 protected:
  Teuchos::ParameterList plist_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
