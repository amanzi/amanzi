/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

The small strain models by Hardin and Drnevich.

*/

#ifndef AMANZI_MECHANICS_SSM_HARDIN_DRNEVICH_HH_
#define AMANZI_MECHANICS_SSM_HARDIN_DRNEVICH_HH_

#include "Teuchos_ParameterList.hpp"

#include "SSM.hh"

namespace Amanzi {
namespace Mechanics {

class SSM_HardinDrnevich : public SSM {
 public:
  SSM_HardinDrnevich(Teuchos::ParameterList& plist)
  {
    Gmax_ = plist.get<double>("maximum shear stress");
    gamma_ref_ = plist.get<double>("reference shear strain");
  }
  virtual ~SSM_HardinDrnevich() {};

  // gamma - shear strain, e - volumetric strain
  virtual double ShearStress(double gamma) { return Gmax_ / (1.0 + gamma / gamma_ref_); }
  virtual double BulkModulus(double e) { return 0.0; }

 private:
  double Gmax_, gamma_ref_;
};

} // namespace Mechanics
} // namespace Amanzi

#endif
