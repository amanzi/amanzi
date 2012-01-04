#ifndef AMANZI_INDICATOR_FACTORY_HH_
#define AMANZI_INDICATOR_FACTORY_HH_

#include <string>
#include <fstream>
#include "Epetra_Comm.h"

#include "indicator.hh"

namespace Amanzi {

class IndicatorFactory {
 public:
  IndicatorFactory() {}
  ~IndicatorFactory() {}
  Indicator* Create(std::string&, const Epetra_Comm&) const;
 private:
  Indicator* create_grid_indicator(int, std::fstream&, const Epetra_Comm&) const;
};

} // namespace Amanzi

#endif // AMANZI_INDICATOR_FACTORY_HH_
