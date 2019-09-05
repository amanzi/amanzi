#include "FunctionPointer.hh"

namespace Amanzi {

FunctionPointer::FunctionPointer(
  double (*f)(const Kokkos::View<double*>&, const Kokkos::View<double*>&), 
  const Kokkos::View<double*> &p)
    : f_(f), np_(0)
{
  if (p.extent(0) > 0) {
    np_ = p.extent(0);
    Kokkos::resize(p_,np_); 
    //p_ = new double[np_];
    for (int i = 0; i < np_; ++i) p_[i] = p[i];
  }
}

FunctionPointer::FunctionPointer(const FunctionPointer& source)
    : f_(source.f_)
{
  if (source.np_) {
    np_ = source.np_;
    Kokkos::resize(p_,np_); 
    //p_ = new double[np_];
    for (int i = 0; i < np_; ++i) p_[i] = source.p_[i];
  }
}

} // namespace Amanzi
