/*
  Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov
*/

#include "FunctionTabularString.hh"
#include "errors.hh"

namespace Amanzi {

FunctionTabularString::FunctionTabularString(const std::vector<double>& x,
                                             const std::vector<std::string>& y)
  : x_(x), y_(y)
{
  CheckArgs_(x, y);
}


void
FunctionTabularString::CheckArgs_(const std::vector<double>& x,
                                  const std::vector<std::string>& y) const
{
  if (x.size() != y.size()) {
    Errors::Message m;
    m << "the number of x and y values differ";
    Exceptions::amanzi_throw(m);
  }
  if (x.size() < 2) {
    Errors::Message m;
    m << "at least two table values must be given";
    Exceptions::amanzi_throw(m);
  }
  for (int j = 1; j < x.size(); ++j) {
    if (x[j] <= x[j - 1]) {
      Errors::Message m;
      m << "x values are not strictly increasing";
      Exceptions::amanzi_throw(m);
    }
  }
}


std::string
FunctionTabularString::operator()(double xv) const
{
  int n = x_.size();
  std::string y;

  if (xv <= x_[0]) {
    y = y_[0];
  } else if (xv > x_[n - 1]) {
    y = y_[n - 1];
  } else {
    // binary search to find interval containing xv
    int j1(0), j2(n - 1);

    while (j2 - j1 > 1) {
      int j = (j1 + j2) / 2;
      // if (xv >= x_[j]) { // right continuous
      if (xv > x_[j]) { // left continuous
        j1 = j;
      } else {
        j2 = j;
      }
    }

    // Now have x_[j1] <= xv < x_[j2], if right continuous
    // or x_[j1] < xv <= x_[j2], if left continuous
    y = y_[j1];
  }

  return y;
}

} // namespace Amanzi
