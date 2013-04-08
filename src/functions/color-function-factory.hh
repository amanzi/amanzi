#ifndef AMANZI_COLOR_FUNCTION_FACTORY_HH_
#define AMANZI_COLOR_FUNCTION_FACTORY_HH_

#include <string>
#include <fstream>
#include "Epetra_Comm.h"

#include "color-function.hh"

namespace Amanzi {

class ColorFunctionFactory {
 public:
  ColorFunctionFactory() {}
  ~ColorFunctionFactory() {}
  ColorFunction* Create(std::string&, const Epetra_Comm&) const;
 private:
  ColorFunction* create_grid_color_function(int, std::fstream&, const Epetra_Comm&) const;
};

} // namespace Amanzi

#endif // AMANZI_COLOR_FUNCTION_FACTORY_HH_
