/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "FunctionPointer.hh"

namespace Amanzi {

FunctionPointer::FunctionPointer(double (*f)(const double*, const double*),
                                 const std::vector<double>& p)
  : f_(f), np_(0), p_(0)
{
  if (p.size() > 0) {
    np_ = p.size();
    p_ = new double[np_];
    for (int i = 0; i < np_; ++i) p_[i] = p[i];
  }
}

FunctionPointer::FunctionPointer(const FunctionPointer& source) : f_(source.f_), np_(0), p_(0)
{
  if (source.p_) {
    np_ = source.np_;
    p_ = new double[np_];
    for (int i = 0; i < np_; ++i) p_[i] = source.p_[i];
  }
}

FunctionPointer::~FunctionPointer()
{
  if (p_) delete[] p_;
}

} // namespace Amanzi
