/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! FunctionDistance: distance from a reference point.
/*!

A distance function calculates distance from reference point :math:`x_0`
using by the following expression:

.. math::
  f(x) = \sum_{j=0}^{n} m_j (x_j - x_{0,j})^2

Note that the first parameter in :math:`x` can be time.
Here is an example of a distance function using isotropic metric:

Example:
.. code-block:: xml

  <ParameterList name="function-squaredistance">
    <Parameter name="x0" type="Array(double)" value="{1.0, 3.0, 0.0}"/>
    <Parameter name="metric" type="Array(double)" value="{1.0, 1.0, 1.0}"/>
  </ParameterList>

*/

#ifndef AMANZI_SQUAREDISTANCE_FUNCTION_HH_
#define AMANZI_SQUAREDISTANCE_FUNCTION_HH_

#include <vector>

#include "Function.hh"

namespace Amanzi {

class FunctionSquareDistance : public Function {
 public:
  FunctionSquareDistance(const std::vector<double>& x0, const std::vector<double>& metric);
  ~FunctionSquareDistance() {}
  std::unique_ptr<Function> Clone() const
  {
    return std::make_unique<FunctionSquareDistance>(*this);
  }
  double operator()(const std::vector<double>& x) const;

 private:
  std::vector<double> x0_, metric_;
};

} // namespace Amanzi

#endif
