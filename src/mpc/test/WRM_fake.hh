/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#ifndef AMANZI_FAKE_MODEL_HH_
#define AMANZI_FAKE_MODEL_HH_

#include "Factory.hh"

#include "WRM.hh"

namespace Amanzi {
namespace Flow {

class WRM_fake : public WRM {
 public:
  explicit WRM_fake(Teuchos::ParameterList& plist)
  {
    alpha = 1.0;
    n = 2.0;
    m = 1.0;
  }
  ~WRM_fake(){};

  // required methods from the base class
  // -- relative permeability formula
  double k_relative(double pc) const { return 1.0; }

  // -- analytic solution was designed for stationary PDE
  double saturation(double pc) const { return -pc; }
  double dSdPc(double pc) const { return -1.0; }

  double capillaryPressure(double saturation) const { return -saturation; }
  double residualSaturation() const { return 0.0; }

  // -- derivative of rel_perm  w.r.t. capillary pressure.
  double dKdPc(double pc) const { return 0.0; }

 private:
  double m, n, alpha;

  static Utils::RegisteredFactory<WRM, WRM_fake> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
