#include "tabular-function.hh"
#include "errors.hh"

namespace Amanzi {

TabularFunction::TabularFunction(const std::vector<double> &x, const std::vector<double> &y)
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
    if (x[j] <= x[j-1]) {
      Errors::Message m;
      m << "x values are not strictly increasing";
      Exceptions::amanzi_throw(m);
    }
  }
  x_ = x;
  y_ = y;
}

double TabularFunction::operator() (const double *x) const
{
  double y;
  int n = x_.size();
  if (*x <= x_[0]) {
    y = y_[0];
  } else if (*x >= x_[n-1]) {
    y = y_[n-1];
  } else {
    // binary search to find interval containing *x
    int j1 = 0, j2 = n-1;
    while (j2 - j1 > 1) {
      int j = (j1 + j2) / 2;
      if (*x > x_[j]) {
        j1 = j;
      } else {
        j2 = j;
      }
    }
    // Linear interpolation between x[j1] and x[j2]
    y = y_[j1] + ((y_[j2]-y_[j1])/(x_[j2]-x_[j1])) * (*x - x_[j1]);
  }
  return y;
}

} // namespace Amanzi
