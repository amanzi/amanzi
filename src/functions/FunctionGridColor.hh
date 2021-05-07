#ifndef AMANZI_GRID_COLOR_FUNCTION_HH_
#define AMANZI_GRID_COLOR_FUNCTION_HH_

#include <vector>

#include "FunctionColor.hh"

namespace Amanzi {

class FunctionGridColor : public FunctionColor {
 public:
  FunctionGridColor(int dim, std::vector<int> &count, std::vector<double> &x0,
                std::vector<double> &dx, std::vector<int> &array)
      : dim_(dim), count_(count), x0_(x0), dx_(dx), array_(array) {}
  ~FunctionGridColor() {}
  std::unique_ptr<FunctionColor> Clone() const {
    return std::make_unique<FunctionGridColor>(*this);
  }
  int operator()(const double* ) const;

private:
  int dim_;
  std::vector<int> count_;
  std::vector<double> x0_;
  std::vector<double> dx_;
  std::vector<int> array_;
};

} // namespace Amanzi

#endif // AMANZI_GRID_COLOR_FUNCTION_HH_
