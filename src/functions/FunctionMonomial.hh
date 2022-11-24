/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! FunctionMonomial: a multivariate monomial function.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

A multi-variable monomial function is given by the following expression:

.. math::
  f(x) = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}

with the constant factor :math:`c`, the reference point :math:`x_0`, and
integer exponents :math:`p_j`. 
Note that the first parameter in :math:`x` can be time.

* `"c`" ``[double]`` c in f = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}
* `"x0`" ``[Array(double)]`` x0 in f = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}
* `"exponents`" ``[Array(int)]`` p in f = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}

Conditions:

.. code-block:: python

  len(x0) == len(exponents)

Here is an example of monomial of degree 6 in three variables:

.. code-block:: xml

  <ParameterList name="function-monomial">
    <Parameter name="c" type="double" value="1.0"/>
    <Parameter name="x0" type="Array(double)" value="{1.0, 3.0, 0.0}"/>
    <Parameter name="exponents" type="Array(int)" value="{2, 3, 1}"/>
  </ParameterList>

*/

#ifndef AMANZI_MONOMIAL_FUNCTION_HH_
#define AMANZI_MONOMIAL_FUNCTION_HH_

#include <vector>

#include "Function.hh"

namespace Amanzi {

class FunctionMonomial : public Function {
 public:
  FunctionMonomial(double c, const std::vector<double>& x0, const std::vector<int>& p);
  ~FunctionMonomial() {}
  std::unique_ptr<Function> Clone() const { return std::make_unique<FunctionMonomial>(*this); }
  double operator()(const std::vector<double>& x) const;

 private:
  double c_;
  std::vector<double> x0_;
  std::vector<int> p_;
};

} // namespace Amanzi

#endif
