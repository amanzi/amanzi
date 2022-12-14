/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1)
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#ifndef AMANZI_VAN_GENUCHTEN_MODEL_HH_
#define AMANZI_VAN_GENUCHTEN_MODEL_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "WRM.hh"

namespace Amanzi {
namespace Flow {

class WRM_vanGenuchten : public WRM {
 public:
  explicit WRM_vanGenuchten(Teuchos::ParameterList& plist);

  explicit WRM_vanGenuchten(double m,
                            double l,
                            double alpha,
                            double sr,
                            std::string& krel_function,
                            double pc0 = 0.0);
  ~WRM_vanGenuchten(){};

  // required methods from the base class
  double k_relative(double pc) const;
  double saturation(double pc) const;
  double dSdPc(double pc) const;
  double capillaryPressure(double saturation) const;
  double residualSaturation() const { return sr_; }
  double dKdPc(double pc) const;

 private:
  void Init_(double m, double l, double alpha, double sr, std::string& krel_function, double pc0);

 private:
  double m_, n_, l_, alpha_; // van Genuchten parameters
  double sr_;                // residual saturation
  int function_;             // relative permeability model
  double tol_;               // defines when cut off derivative which tends to go to infinity

  double pc0_;                  // regularization threshold (usually 0 to 500 Pa)
  double a_, b_, factor_dSdPc_; // frequently used constant

  static Utils::RegisteredFactory<WRM, WRM_vanGenuchten> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
