/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

The list of multiscale models is the place for various subscale models that coul 
be mixed and matched.
Its ordered by materials and includes parameters for the assigned multiscale model
This list is optional.

* `"multiscale model`" [string] is the model name. Available option is `"dual porosity`"
  and `"generalized dual porosity`".

* `"regions`" [Array(string)] is the list of regions where this model should be applied.

* `"xxx parameters`" [sublist] provides parameters for the model specified by variable `"multiscale model`".

*/

#ifndef MULTISCALE_TRANSPORT_POROSITY_HH_
#define MULTISCALE_TRANSPORT_POROSITY_HH_

#include <string>
#include <vector>

#include "DenseVector.hh"

namespace Amanzi {
namespace Transport {

class MultiscaleTransportPorosity {
 public:
  virtual ~MultiscaleTransportPorosity(){};

  // Compute solute flux: icomp - component id, phi - matrix porosity,
  // tcc_m_aux - vector of concentration values in secondary nodes,
  // wfm[0|1] - fracture water content at initial and final time moments,
  // wcm[0|1] - matrix water content at initial and final time moments
  virtual double ComputeSoluteFlux(double flux_liquid,
                                   double& tcc_f,
                                   WhetStone::DenseVector& tcc_m,
                                   int icomp,
                                   double dt,
                                   double wcf0,
                                   double wcf1,
                                   double wcm0,
                                   double wcm1,
                                   double phi) = 0;

  // Number of matrix nodes
  virtual int NumberMatrixNodes() = 0;
};

} // namespace Transport
} // namespace Amanzi

#endif
