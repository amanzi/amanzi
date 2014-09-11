#ifndef AMANZI_FUNCTION_HH_
#define AMANZI_FUNCTION_HH_

#include <vector>

namespace Amanzi {

class Function {
 public:
  virtual ~Function() {}
  virtual Function* Clone() const = 0;
  virtual double operator()(const std::vector<double>& ) const = 0;
};

} // namespace Amanzi

#endif // AMANZI_FUNCTION_HH_
