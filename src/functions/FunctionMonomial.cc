/*
  Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <math.h>

#include "errors.hh"
#include "FunctionMonomial.hh"

namespace Amanzi {

FunctionMonomial::FunctionMonomial(double c,
                                   const std::vector<double>& x0,
                                   const std::vector<int>& p)
{
  if (x0.size() != p.size()) {
    Errors::Message m;
    m << "Mismatch of multi-index and reference point dimensions.";
    Exceptions::amanzi_throw(m);
  }
  c_ = c;
  x0_ = x0;
  p_ = p;
}


double
FunctionMonomial::operator()(const std::vector<double>& x) const
{
  double y = c_;
  if (x.size() < x0_.size()) {
    Errors::Message m;
    m << "FunctionMonomial expects higher-dimensional argument.";
    Exceptions::amanzi_throw(m);
  }
  for (int j = 0; j < x0_.size(); ++j) y *= pow(x[j] - x0_[j], p_[j]);
  return y;
}

} // namespace Amanzi
