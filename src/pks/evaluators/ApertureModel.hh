/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

Aperture-stress relations like Batron-Bandis' equation are used to estimate
changes in fracture aperture due to change in internal fluid/gas pressure.

*/

#ifndef AMANZI_EVALUATORS_APERTURE_MODEL_HH_
#define AMANZI_EVALUATORS_APERTURE_MODEL_HH_

#include <string>

namespace Amanzi {
namespace Evaluators {

class ApertureModel {
 public:
  virtual ~ApertureModel() {};
  virtual double OverburdenPressure() = 0;
  virtual double Aperture(double pressure, double pov) = 0;
  virtual double dAperturedPressure(double pressure, double pov) = 0;
};

} // namespace Evaluators
} // namespace Amanzi

#endif
