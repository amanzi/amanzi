#include "pointer-function.hh"

namespace Amanzi {

PointerFunction::PointerFunction(double (*f)(const double*, const double*), const std::vector<double> &p)
    : f_(f), np_(0), p_(0)
{
  if (p.size() > 0) {
    np_ = p.size();
    p_ = new double[np_];
    for (int i = 0; i < np_; ++i) p_[i] = p[i];
  }
}

PointerFunction::PointerFunction(const PointerFunction& source)
    : f_(source.f_), np_(0), p_(0)
{
  if (source.p_) {
    np_ = source.np_;
    p_ = new double[np_];
    for (int i = 0; i < np_; ++i) p_[i] = source.p_[i];
  }
}

PointerFunction::~PointerFunction()
{
  if (p_) delete [] p_;
}

} // namespace Amanzi
