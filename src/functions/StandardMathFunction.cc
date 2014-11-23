#include <math.h>

#include "StandardMathFunction.hh"
#include "errors.hh"

namespace Amanzi {

StandardMathFunction::StandardMathFunction(std::string op,
        double amplitude, double parameter, double shift) :
    op_(op), amplitude_(amplitude), parameter_(parameter), shift_(shift) {
  if (!((op_ == "cos") ||
        (op_ == "sin") ||
        (op_ == "tan") ||
        (op_ == "acos") ||
        (op_ == "asin") ||
        (op_ == "atan") ||
        (op_ == "cosh") ||
        (op_ == "sinh") ||
        (op_ == "tanh") ||
        (op_ == "exp") ||
        (op_ == "log") ||
        (op_ == "log10") ||
        (op_ == "sqrt") ||
        (op_ == "ceil") ||
        (op_ == "fabs") ||
        (op_ == "floor") ||
        (op_ == "mod") ||
        (op_ == "pow"))) {
    std::stringstream m;
    m << "Invalid or unknown standard math function " << op_;
    Errors::Message message(m.str());
    Exceptions::amanzi_throw(message);
  }
}

double StandardMathFunction::operator()(const std::vector<double>& x) const
{
  double x0 = x[0] - shift_;
  if (op_ == "cos") {
    return amplitude_ * cos(x0);
  } else if (op_ == "sin") {
    return amplitude_ * sin(x0);
  } else if (op_ == "tan") {
    return amplitude_ * tan(x0);
  } else if (op_ == "acos") {
    return amplitude_ * acos(x0);
  } else if (op_ == "asin") {
    return amplitude_ * asin(x0);
  } else if (op_ == "atan") {
    return amplitude_ * atan(x0);
  } else if (op_ == "cosh") {
    return amplitude_ * cosh(x0);
  } else if (op_ == "sinh") {
    return amplitude_ * sinh(x0);
  } else if (op_ == "tanh") {
    return amplitude_ * tanh(x0);
  } else if (op_ == "exp") {
    return amplitude_ * exp(x0);
  } else if (op_ == "log") {
    if (x0 <= 0) InvalidDomainError_(x[0]);
    return amplitude_ * log(x0);
  } else if (op_ == "log10") {
    if (x0 <= 0) InvalidDomainError_(x[0]);
    return amplitude_ * log10(x0);
  } else if (op_ == "sqrt") {
    if (x0 < 0) InvalidDomainError_(x[0]);
    return amplitude_ * sqrt(x0);
  } else if (op_ == "ceil") {
    return amplitude_ * ceil(x0);
  } else if (op_ == "fabs") {
    return amplitude_ * fabs(x0);
  } else if (op_ == "floor") {
    return amplitude_ * floor(x0);
  } else if (op_ == "pow") {
    return amplitude_ * pow(x0, parameter_);
  } else if (op_ == "mod") {
    return fmod(x0, parameter_);
  } else {
    std::stringstream m;
    m << "Invalid or unknown standard math function " << op_;
    Errors::Message message(m.str());
    Exceptions::amanzi_throw(message);
  }
  return 0.0;
}

void StandardMathFunction::InvalidDomainError_(double x) const {
  std::stringstream m;
  m << "Value " << x << " is not in the domain of operator " << op_ << ".";
  Errors::Message message(m.str());
  Exceptions::amanzi_throw(message);
}

} // namespace Amanzi
