#ifndef AMANZI_COLOR_FUNCTION_FACTORY_HH_
#define AMANZI_COLOR_FUNCTION_FACTORY_HH_

#include <string>
#include <fstream>

#include "AmanziTypes.hh"
#include "FunctionColor.hh"

namespace Amanzi {

class FunctionColorFactory {
 public:
  FunctionColorFactory() {}
  ~FunctionColorFactory() {}
  std::unique_ptr<FunctionColor> Create(std::string&, const Comm_type&) const;
 private:
  std::unique_ptr<FunctionColor> create_grid_color_function(int, std::fstream&, const Comm_type&) const;
};

} // namespace Amanzi

#endif // AMANZI_COLOR_FUNCTION_FACTORY_HH_
