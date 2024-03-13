/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Functions

*/

#ifndef AMANZI_TABULAR_STRING_FUNCTION_HH_
#define AMANZI_TABULAR_STRING_FUNCTION_HH_

#include <string>
#include <vector>

#include "Function.hh"

namespace Amanzi {

class FunctionTabularString {
 public:
  FunctionTabularString(const std::vector<double>& x, const std::vector<std::string>& y);
  ~FunctionTabularString(){};

  std::string operator()(double xv) const;

 private:
  std::vector<double> x_;
  std::vector<std::string> y_;

 private:
  void CheckArgs_(const std::vector<double>& x, const std::vector<std::string>& y) const;
};

} // namespace Amanzi

#endif
