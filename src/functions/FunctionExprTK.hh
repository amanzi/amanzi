/*
  Functions

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

/*!

This function parses a string expression. The function has min(N, D + 1)
arguments t, x, y, and z. The argument t is required. D is the space dimension,
and N is the user specified number of arguments which could be less than D + 1.

Example:

.. code-block:: xml

  <ParameterList name="function-exprtk">
    <Parameter name="number of arguments" type="int" value="3"/>
    <Parameter name="formula" type="string" value="t + x + 2 * y"/>
  </ParameterList>
*/

#ifndef AMANZI_EXPRTK_FUNCTION_HH_
#define AMANZI_EXPRTK_FUNCTION_HH_

#include <memory>
#include <vector>

#include "ExprTK.hh"

#include "Function.hh"

namespace Amanzi {

class FunctionExprTK : public Function {
 public:
  FunctionExprTK(int n, const std::string& formula);
  FunctionExprTK* Clone() const { return new FunctionExprTK(*this); }
  double operator()(const std::vector<double>& x) const;

 private:
  std::shared_ptr<Utils::ExprTK> exprtk_;
};

} // namespace Amanzi

#endif
