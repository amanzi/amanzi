#ifndef AMANZI_INDICATOR_FACTORY_HH_
#define AMANZI_INDICATOR_FACTORY_HH_

#include <string>
#include <fstream>
#include "Epetra_Comm.h"

namespace Amanzi {

class Indicator; // forward declaration

class IndicatorFactory {
 public:
  IndicatorFactory() {}
  ~IndicatorFactory() {}
  Indicator* Create(std::string&, Epetra_Comm&) const;
 private:
  Indicator* create_grid_indicator(int, std::fstream&, Epetra_Comm&) const;
};

} // namespace Amanzi

#endif // AMANZI_INDICATOR_FACTORY_HH_
