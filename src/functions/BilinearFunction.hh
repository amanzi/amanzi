#ifndef AMANZI_BILINEAR_FUNCTION_HH_
#define AMANZI_BILINEAR_FUNCTION_HH_

#include <vector>

#include "Function.hh"
#include "Epetra_SerialDenseMatrix.h"

namespace Amanzi {

class BilinearFunction : public Function {
 public:
  BilinearFunction(const std::vector<double> &x, const std::vector<double> &y,
                   const Epetra_SerialDenseMatrix &v, const int xi, const int yi);
  ~BilinearFunction() {};
  BilinearFunction* Clone() const { return new BilinearFunction(*this); }
  double operator()(const std::vector<double>& x) const;

 private:
  std::vector<double> x_, y_;
  Epetra_SerialDenseMatrix v_;
  int xi_, yi_;
  
 private: // helper functions
  void check_args(const std::vector<double>&, const std::vector<double>&, const Epetra_SerialDenseMatrix&) const;
};

} // namespace Amanzi

#endif  // AMANZI_BILINEAR_FUNCTION_HH_
