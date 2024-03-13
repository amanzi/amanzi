/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef AMANZI_COLOR_FUNCTION_HH_
#define AMANZI_COLOR_FUNCTION_HH_

#include <memory>

namespace Amanzi {

class FunctionColor {
 public:
  virtual ~FunctionColor() {}
  virtual std::unique_ptr<FunctionColor> Clone() const = 0;
  virtual int operator()(const double*) const = 0;

  virtual int getDimension() const = 0;
};

} // namespace Amanzi

#endif //  AMANZI_COLOR_FUNCTION_HH_
