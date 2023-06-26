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

#ifndef AMANZI_LINEAR_MODEL_HH_
#define AMANZI_LINEAR_MODEL_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "WRM.hh"

namespace Amanzi {
namespace Flow {

class WRM_linear : public WRM {
 public:
  explicit WRM_linear(Teuchos::ParameterList& plist);

  explicit WRM_linear(double pc0);
  ~WRM_linear(){};

  // required methods from the base class
  double k_relative(double pc) const;
  double saturation(double pc) const;
  double dSdPc(double pc) const;
  double capillaryPressure(double saturation) const;
  double residualSaturation() const { return 0.0; }
  double dKdPc(double pc) const;

 private:
  double pc0_, sr_;

  static Utils::RegisteredFactory<WRM, WRM_linear> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
