#include "FunctionTabular.hh"
#include "errors.hh"

namespace Amanzi {

FunctionTabular::FunctionTabular(const Kokkos::View<double*>& x, const Kokkos::View<double*>& y,
                                 const int xi)
  : x_(x), y_(y), xi_(xi)
{
  Kokkos::resize(form_,x.extent(0)-1);
  for(int i = 0 ; i < form_.extent(0); ++i){
    form_(i) = LINEAR; 
  } 
  //form_.assign(x.size() - 1, LINEAR);
  check_args(x, y, form_);
}

FunctionTabular::FunctionTabular(
    const Kokkos::View<double*>& x, const Kokkos::View<double*>& y,
    const int xi, const Kokkos::View<Form*>& form) : x_(x), y_(y), xi_(xi), form_(form)
{
  check_args(x, y, form);
}

FunctionTabular::FunctionTabular(
    const Kokkos::View<double*>& x, const Kokkos::View<double*>& y,
    const int xi, const Kokkos::View<Form*>& form, const std::vector<Function*>& func) 
  : x_(x), y_(y), xi_(xi), form_(form), func_(func)
{
  check_args(x, y, form);
}

void FunctionTabular::check_args(const Kokkos::View<double*>& x, const Kokkos::View<double*>& y,
                                 const Kokkos::View<Form*>& form) const
{
  if (x.extent(0) != y.extent(0)) {
    Errors::Message m;
    m << "the number of x and y values differ";
    Exceptions::amanzi_throw(m);
  }
  if (x.extent(0) < 2) {
    Errors::Message m;
    m << "at least two table values must be given";
    Exceptions::amanzi_throw(m);
  }
  for (int j = 1; j < x.extent(0); ++j) {
    if (x[j] <= x[j-1]) {
      Errors::Message m;
      m << "x values are not strictly increasing";
      Exceptions::amanzi_throw(m);
    }
  }
  if (form.extent(0) != x.extent(0)-1) {
    Errors::Message m;
    m << "incorrect number of form values specified";
    Exceptions::amanzi_throw(m);
  }
}

double FunctionTabular::operator()(const Kokkos::View<double*>& x) const
{
  double y;
  double xv = x[xi_];
  int n = x_.extent(0);
  if (xv <= x_[0]) {
    y = y_[0];
  } else if (xv > x_[n-1]) {
    y = y_[n-1];
  } else {
    // binary search to find interval containing xv
    int j1 = 0, j2 = n-1;
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
      y = y_[j1] + ((y_[j2]-y_[j1])/(x_[j2]-x_[j1])) * (xv - x_[j1]);
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
