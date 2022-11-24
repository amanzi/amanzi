#include <algorithm>

#include "FunctionPolynomial.hh"
#include "errors.hh"

namespace Amanzi {

FunctionPolynomial::FunctionPolynomial(const std::vector<double>& c,
                                       const std::vector<int>& p,
                                       double x0)
{
  if (c.size() < 1) {
    Errors::Message m;
    m << "at least one value requred for the coefficient vector";
    Exceptions::amanzi_throw(m);
  }
  if (p.size() != c.size()) {
    Errors::Message m;
    m << "the number of values for the coefficient and exponent vectors differ";
    Exceptions::amanzi_throw(m);
  }
  // Minimum and maximum powers.
  pmin_ = std::min(0, *(std::min_element(p.begin(), p.end())));
  pmax_ = std::max(0, *(std::max_element(p.begin(), p.end())));
  int n = pmax_ - pmin_ + 1;
  c_.resize(n);
  c_.assign(n, 0.0);
  for (int j = 0; j < c.size(); ++j) c_[p[j] - pmin_] += c[j];
  x0_ = x0;
}

double
FunctionPolynomial::operator()(const std::vector<double>& x) const
{
  // Polynomial terms with non-negative exponents
  double y = c_[pmax_ - pmin_];
  if (pmax_ > 0) {
    double z = x[0] - x0_;
    for (int j = pmax_; j > 0; --j) y = c_[j - 1 - pmin_] + z * y;
  }
  // Polynomial terms with negative exponents.
  if (pmin_ < 0) {
    double w = c_[0];
    double z = 1.0 / (x[0] - x0_);
    for (int j = pmin_; j < -1; ++j) w = c_[j + 1 - pmin_] + z * w;
    y += z * w;
  }
  return y;
}

} // namespace Amanzi
