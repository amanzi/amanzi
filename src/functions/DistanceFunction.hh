/*
  Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_DISTANCE_FUNCTION_HH_
#define AMANZI_DISTANCE_FUNCTION_HH_

#include <vector>

#include "Function.hh"

namespace Amanzi {

class DistanceFunction : public Function {
 public:
  DistanceFunction(const std::vector<double> &x0);
  ~DistanceFunction() {}
  DistanceFunction* Clone() const { return new DistanceFunction(*this); }
  double operator()(const std::vector<double>& x) const;

 private:
  std::vector<double> x0_;
};

} // namespace Amanzi

#endif
