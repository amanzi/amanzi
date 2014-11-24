#ifndef AMANZI_STANDARD_MATH_FUNCTION_HH_
#define AMANZI_STANDARD_MATH_FUNCTION_HH_

#include <vector>
#include <string>

#include "Function.hh"

namespace Amanzi {

class StandardMathFunction : public Function {

public:
  StandardMathFunction(std::string op, double amplitude, double parameter, double shift);
  ~StandardMathFunction() {}
  StandardMathFunction* Clone() const { return new StandardMathFunction(*this); }
  double operator()(const std::vector<double>& x) const;

private:
  void InvalidDomainError_(double x) const;

private:
  double parameter_;
  double amplitude_;
  double shift_;
  std::string op_;
};

} // namespace Amanzi

#endif  // AMANZI_STANDARD_MATH_FUNCTION_HH_
