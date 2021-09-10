/*
  Flow PK
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
  We use this class to test convergence of discretization
  schemes. It employs a simple model for relative permeability,
  k_rel = 1 / (1 + pc^2).
*/

#ifndef AMANZI_FAKE_MODEL_HH_
#define AMANZI_FAKE_MODEL_HH_

#include "Factory.hh"

#include "WRM.hh"

namespace Amanzi {
namespace Flow {

class WRM_fake : public WRM {
 public:
  explicit WRM_fake(Teuchos::ParameterList& plist) {
    alpha = 1.0;
    n = 2.0;
    m = 1.0;
  }
  ~WRM_fake() {};
  
  // required methods from the base class
  // -- relative permeability formula
  double k_relative(double pc) const {
    if (pc < 0.0)
      return 1.0 / (1.0 + pc * pc);
    else
      return 1.0;
  }

  // -- analytic solution was designed for stationary PDE
  double saturation(double pc) const { return -pc; }
  double dSdPc(double pc) const { return -1.0; }

  double capillaryPressure(double saturation) const { return -saturation; }
  double residualSaturation() const { return 0.0; }

  // -- derivative of rel_perm  w.r.t. capillary pressure.
  double dKdPc(double pc) const {
    if (pc < 0.0) {
      double tmp = 1.0 + pc * pc;
      return -2 * pc / (tmp * tmp);
    } else {
      return 0.0;
    }
  }

 private:
  double m, n, alpha;

  static Utils::RegisteredFactory<WRM, WRM_fake> factory_;
};

}  // namespace Flow
}  // namespace Amanzi
 
#endif
