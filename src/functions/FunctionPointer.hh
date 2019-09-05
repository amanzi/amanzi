#ifndef AMANZI_POINTER_FUNCTION_HH_
#define AMANZI_POINTER_FUNCTION_HH_

#include <vector>

#include "Function.hh"

namespace Amanzi {

class FunctionPointer : public Function {
 public:
  FunctionPointer(double (*f)(const Kokkos::View<double*>&, const Kokkos::View<double*>&)) : f_(f), np_(0) {}
  FunctionPointer(double (*f)(const Kokkos::View<double*>&, const Kokkos::View<double*>&), const Kokkos::View<double*>&);
  FunctionPointer(const FunctionPointer&);
  ~FunctionPointer(){};
  FunctionPointer* Clone() const { return new FunctionPointer(*this); }
  double operator()(const Kokkos::View<double*>& x) const { return (*f_)(x, p_); }

  KOKKOS_INLINE_FUNCTION double apply_gpu(const Kokkos::View<double*>& x) const {assert(false); return 0.0;}
  
  void apply(const Kokkos::View<double*>& in, Kokkos::View<double*>& out){}

 private:
  double (*f_)(const Kokkos::View<double*>&, const Kokkos::View<double*>&);
  int np_;
  Kokkos::View<double*> p_;
};

} // namespace Amanzi

#endif // AMANZI_POINTER_FUNCTION_HH_
