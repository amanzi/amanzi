/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

Relations like the Kozeny–Carman equation are used to estimate changes in permeability
due to changes in porosity, e.g. due to mineral precipitation or biofilm accumulation.
The change in permeability is calculated by relating the current permeability, :math:`K`,
based on the current porosity, :math:`\phi`, to the initial or reference permeability
and porosity:

.. math::
  K = K_0 \frac{f(\phi)}{f(\phi_0)}

*/

#ifndef AMANZI_FLOW_PERMEABILITY_HH_
#define AMANZI_FLOW_PERMEABILITY_HH_

#include <string>

namespace Amanzi {
namespace Flow {

class Permeability {
 public:
  virtual ~Permeability(){};
  virtual double Factor(double phi) = 0;
  virtual double dFactordPorosity(double phi) = 0;
};

} // namespace Flow
} // namespace Amanzi

#endif
