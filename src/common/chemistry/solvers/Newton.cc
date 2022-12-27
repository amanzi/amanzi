/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry

*/

#include "Newton.hh"

namespace Amanzi {
namespace AmanziChemistry {

Newton::Newton(const int n)
{
  size(n);
  x_.resize(n);
  r_.resize(n);
  indices_.resize(n);
  vv_.resize(n);
}

void
Newton::solve()
{
  std::cout << "Solved!\n";
}

} // namespace AmanziChemistry
} // namespace Amanzi
