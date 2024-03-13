/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

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
  std::unique_ptr<FunctionColor>
  create_grid_color_function(int, std::fstream&, const Comm_type&) const;
};

} // namespace Amanzi

#endif // AMANZI_COLOR_FUNCTION_FACTORY_HH_
