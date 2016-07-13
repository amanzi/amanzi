/*
  Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "DistanceFunction.hh"
#include "errors.hh"

namespace Amanzi {

DistanceFunction::DistanceFunction(const std::vector<double>& x0)
{
  x0_ = x0;
}


double DistanceFunction::operator()(const std::vector<double>& x) const
{
  double tmp, y(0.0);
  if (x.size() < x0_.size()) {
    Errors::Message m;
    m << "DistanceFunction expects higher-dimensional argument.";
    Exceptions::amanzi_throw(m);
  }    
  for (int j = 0; j < x0_.size(); ++j) {
    tmp = x[j] - x0_[j];
    y += tmp * tmp;
  }
  return y;
}

}  // namespace Amanzi
