/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//! Function: base class for all functions of space and time.
/*!
Analytic, algabraic functions of space and time are used for a variety of
purposes, including boundary conditions, initial conditions, and independent
variables.

For initial conditions, functions are prescribed of space only, i.e.

:math:`u = f(x,y,z)`

For boundary conditions and independent variables, functions are also a
function of time:

:math:`u = f(t,x,y,z)`

Note, this does not follow the `"typed`" format for legacy reasons.

*/

#ifndef AMANZI_FUNCTION_HH_
#define AMANZI_FUNCTION_HH_

#include <vector>
#include <memory>

namespace Amanzi {

enum class Form_kind { LINEAR, CONSTANT, FUNCTION };

class Function {
 public:
  virtual ~Function() {}
  virtual std::unique_ptr<Function> Clone() const = 0;
  virtual double operator()(const std::vector<double>&) const = 0;
};

} // namespace Amanzi

#endif // AMANZI_FUNCTION_HH_
