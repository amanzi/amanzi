/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef AMANZI_POINTER_FUNCTION_HH_
#define AMANZI_POINTER_FUNCTION_HH_

#include <vector>

#include "Function.hh"

namespace Amanzi {

class FunctionPointer : public Function {
 public:
  FunctionPointer(double (*f)(const double*, const double*)) : f_(f), np_(0), p_(0) {}
  FunctionPointer(double (*f)(const double*, const double*), const std::vector<double>&);
  FunctionPointer(const FunctionPointer&);
  ~FunctionPointer();
  std::unique_ptr<Function> Clone() const { return std::make_unique<FunctionPointer>(*this); }
  double operator()(const std::vector<double>& x) const { return (*f_)(&x[0], p_); }

 private:
  double (*f_)(const double*, const double*);
  int np_;
  double* p_;
};

} // namespace Amanzi

#endif // AMANZI_POINTER_FUNCTION_HH_
