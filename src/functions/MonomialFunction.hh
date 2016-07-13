/*
  Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MONOMIAL_FUNCTION_HH_
#define AMANZI_MONOMIAL_FUNCTION_HH_

#include <vector>

#include "Function.hh"

namespace Amanzi {

class MonomialFunction : public Function {
 public:
  MonomialFunction(double c, const std::vector<double> &x0, const std::vector<int> &p);
  ~MonomialFunction() {}
  MonomialFunction* Clone() const { return new MonomialFunction(*this); }
  double operator()(const std::vector<double>& x) const;

 private:
  double c_;
  std::vector<double> x0_;
  std::vector<int> p_;
};

} // namespace Amanzi

#endif
