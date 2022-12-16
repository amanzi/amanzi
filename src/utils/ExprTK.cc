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

#include "ExprTK.hh"

namespace Amanzi {
namespace Utils {

bool
ExprTK::Initialize(int n, const std::string& formula)
{
  n_ = n;
  symbol_table_.add_variable("t", t);
  if (n_ > 1) symbol_table_.add_variable("x", x);
  if (n_ > 2) symbol_table_.add_variable("y", y);
  if (n_ > 3) symbol_table_.add_variable("z", z);

  t = x = y = z = 0.0;

  expression_.register_symbol_table(symbol_table_);

  if (!parser_.compile(formula, expression_)) return false;
  return true;
}


double
ExprTK::operator()(const std::vector<double>& txyz)
{
  t = txyz[0];
  if (n_ > 1) x = txyz[1];
  if (n_ > 2) y = txyz[2];
  if (n_ > 3) z = txyz[3];
  return expression_.value();
}

} // namespace Utils
} // namespace Amanzi
