#ifndef AMANZI_COLOR_FUNCTION_HH_
#define AMANZI_COLOR_FUNCTION_HH_

#include <memory>
#include "UniqueHelpers.hh"

namespace Amanzi {

class FunctionColor {
 public:
  virtual ~FunctionColor() {}
  virtual std::unique_ptr<FunctionColor> Clone() const = 0;
  virtual int operator()(const double* ) const = 0;
};

} // namespace Amanzi

#endif //  AMANZI_COLOR_FUNCTION_HH_
