#ifndef AMANZI_COLOR_FUNCTION_HH_
#define AMANZI_COLOR_FUNCTION_HH_

namespace Amanzi {

class ColorFunction {
  public:
    virtual ~ColorFunction() {}
    virtual ColorFunction* Clone() const = 0;
    virtual int operator() (const double*) const = 0;
};

} // namespace Amanzi

#endif //  AMANZI_COLOR_FUNCTION_HH_
