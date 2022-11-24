/*
  Copyright 2010-202X held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

#include "FunctionTabular.hh"
#include "errors.hh"

namespace Amanzi {

FunctionTabular::FunctionTabular(const std::vector<double>& x,
                                 const std::vector<double>& y,
                                 const int xi)
  : x_(x), y_(y), xi_(xi)
{
  form_.assign(x.size() - 1, LINEAR);
  check_args(x, y, form_);
}

FunctionTabular::FunctionTabular(const std::vector<double>& x,
                                 const std::vector<double>& y,
                                 const int xi,
                                 const std::vector<Form>& form)
  : x_(x), y_(y), xi_(xi), form_(form)
{
  check_args(x, y, form);
}

FunctionTabular::FunctionTabular(const std::vector<double>& x,
                                 const std::vector<double>& y,
                                 const int xi,
                                 const std::vector<Form>& form,
                                 std::vector<std::unique_ptr<Function>> func)
  : x_(x), y_(y), xi_(xi), form_(form), func_(std::move(func))
{
  check_args(x, y, form);
}


FunctionTabular::FunctionTabular(const FunctionTabular& other)
  : x_(other.x_), y_(other.y_), xi_(other.xi_), form_(other.form_), func_()
{
  for (const auto& f : other.func_) { func_.emplace_back(f->Clone()); }
}


void
FunctionTabular::check_args(const std::vector<double>& x,
                            const std::vector<double>& y,
                            const std::vector<Form>& form) const
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
  if (form.size() != x.size() - 1) {
    Errors::Message m;
    m << "incorrect number of form values specified";
    Exceptions::amanzi_throw(m);
  }
}

double
FunctionTabular::operator()(const std::vector<double>& x) const
{
  double y;
  double xv = x[xi_];
  int n = x_.size();
  if (xv <= x_[0]) {
    y = y_[0];
  } else if (xv > x_[n - 1]) {
    y = y_[n - 1];
  } else {
    // binary search to find interval containing xv
    int j1 = 0, j2 = n - 1;
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
    switch (form_[j1]) {
    case LINEAR:
      // Linear interpolation between x[j1] and x[j2]
      y = y_[j1] + ((y_[j2] - y_[j1]) / (x_[j2] - x_[j1])) * (xv - x_[j1]);
      break;
    case CONSTANT:
      y = y_[j1];
      break;
    case FUNCTION:
      y = (*func_[j1])(x);
    }
  }
  return y;
}

} // namespace Amanzi
