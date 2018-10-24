#ifndef AMANZI_COLOR_FUNCTION_FACTORY_HH_
#define AMANZI_COLOR_FUNCTION_FACTORY_HH_

#include <string>
#include <fstream>
#include "Epetra_Comm.h"

#include "ColorFunction.hh"

namespace Amanzi {

class ColorFunctionFactory {
 public:
  ColorFunctionFactory() {}
  ~ColorFunctionFactory() {}
  ColorFunction* Create(std::string&, Comm_ptr_type) const;
 private:
  ColorFunction* create_grid_color_function(int, std::fstream&, Comm_ptr_type) const;
};

} // namespace Amanzi

#endif // AMANZI_COLOR_FUNCTION_FACTORY_HH_
