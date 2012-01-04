#ifndef AMANZI_INDICATOR_HH_
#define AMANZI_INDICATOR_HH_

namespace Amanzi {

class Indicator {
  public:
    virtual ~Indicator() {}
    virtual Indicator* Clone() const = 0;
    virtual int operator() (const double*) const = 0;
};

} // namespace Amanzi

#endif //  AMANZI_INDICATOR_HH_
