/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  We use this class for saturated flow, k_rel = 1.
*/

#ifndef AMANZI_SATURATED_MODEL_HH_
#define AMANZI_SATURATED_MODEL_HH_

#include "Factory.hh"

#include "WRM.hh"

namespace Amanzi {
namespace Flow {

class WRM_saturated : public WRM {
 public:
  explicit WRM_saturated(Teuchos::ParameterList& plist){};
  ~WRM_saturated(){};

  // required methods from the base class
  double k_relative(double pc) const { return 1.0; }
  double saturation(double pc) const { return 1.0; }
  double dSdPc(double pc) const { return 0.0; }
  double capillaryPressure(double saturation) const { return 0.0; }
  double residualSaturation() const { return 0.0; }
  double dKdPc(double pc) const { return 0.0; }

 private:
  static Utils::RegisteredFactory<WRM, WRM_saturated> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
