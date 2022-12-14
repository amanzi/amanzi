/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#ifndef AMANZI_WATER_RETENTION_MODEL_HH_
#define AMANZI_WATER_RETENTION_MODEL_HH_

#include <string>

namespace Amanzi {
namespace Flow {

class WRM {
 public:
  virtual ~WRM(){};

  virtual double k_relative(double pc) const = 0;
  virtual double saturation(double pc) const = 0;
  virtual double
  dSdPc(double pc) const = 0; // derivative of saturation w.r.t. to capillary pressure
  virtual double capillaryPressure(double s) const = 0;
  virtual double residualSaturation() const = 0;
  virtual double dKdPc(double pc) const = 0;
};

typedef double (WRM::*KRelFn)(double pc) const;

} // namespace Flow
} // namespace Amanzi

#endif
