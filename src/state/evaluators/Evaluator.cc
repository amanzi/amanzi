/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  State

  Implementation of base class Evaluator helpers.
*/

#include "Evaluator.hh"

namespace Amanzi {
std::ostream&
operator<<(std::ostream& os, const Evaluator& self)
{
  return self.writeInfo(os);
}


} // namespace Amanzi
