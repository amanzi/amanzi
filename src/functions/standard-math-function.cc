#include <math.h>
#include "standard-math-function.hh"
#include "errors.hh"

namespace Amanzi {

StandardMathFunction::StandardMathFunction(std::string op,
        double amplitude, double parameter) :
    op_(op), amplitude_(amplitude), parameter_(parameter) {
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
        (op_ == "pow"))) {
    std::stringstream m;
    m << "Invalid or unknown standard math function " << op_;
    Errors::Message message(m.str());
    Exceptions::amanzi_throw(message);
  }
}

double StandardMathFunction::operator() (const double *x) const
{
  if (op_ == "cos") {
    return amplitude_ * cos(x[0]);
  } else if (op_ == "sin") {
    return amplitude_ * sin(x[0]);
  } else if (op_ == "tan") {
    return amplitude_ * tan(x[0]);
  } else if (op_ == "acos") {
    return amplitude_ * acos(x[0]);
  } else if (op_ == "asin") {
    return amplitude_ * asin(x[0]);
  } else if (op_ == "atan") {
    return amplitude_ * atan(x[0]);
  } else if (op_ == "cosh") {
    return amplitude_ * cosh(x[0]);
  } else if (op_ == "sinh") {
    return amplitude_ * sinh(x[0]);
  } else if (op_ == "tanh") {
    return amplitude_ * tanh(x[0]);
  } else if (op_ == "exp") {
    return amplitude_ * exp(x[0]);
  } else if (op_ == "log") {
    return amplitude_ * log(x[0]);
  } else if (op_ == "log10") {
    return amplitude_ * log10(x[0]);
  } else if (op_ == "sqrt") {
    if (x[0] < 0) InvalidDomainError_(x[0]);
    return amplitude_ * sqrt(x[0]);
  } else if (op_ == "ceil") {
    return amplitude_ * ceil(x[0]);
  } else if (op_ == "fabs") {
    return amplitude_ * fabs(x[0]);
  } else if (op_ == "floor") {
    return amplitude_ * floor(x[0]);
  } else if (op_ == "pow") {
    return amplitude_ * pow(x[0], parameter_);
  } else {
    std::stringstream m;
    m << "Invalid or unknown standard math function " << op_;
    Errors::Message message(m.str());
    Exceptions::amanzi_throw(message);
  }
}

void StandardMathFunction::InvalidDomainError_(double x) const {
  std::stringstream m;
  m << "Value " << x << " is not in the domain of operator " << op_ << ".";
  Errors::Message message(m.str());
  Exceptions::amanzi_throw(message);
}

} // namespace Amanzi
