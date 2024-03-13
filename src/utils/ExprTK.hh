/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella (raovgarimella@lanl.gov)
*/

/*
  Utils

*/

#ifndef AMANZI_EXPRTK_HH_
#define AMANZI_EXPRTK_HH_

#include "exprtk.hpp"

namespace Amanzi {
namespace Utils {

class ExprTK {
 public:
  ExprTK() : n_(0){};
  bool Initialize(int n, const std::string& formula);

  double operator()(const std::vector<double>& txyz);

 private:
  int n_;
  double t, x, y, z;
  exprtk::symbol_table<double> symbol_table_;
  exprtk::expression<double> expression_;
  exprtk::parser<double> parser_;
};

} // namespace Utils
} // namespace Amanzi

#endif
