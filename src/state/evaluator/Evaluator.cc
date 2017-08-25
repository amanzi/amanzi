/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

  Arcos

  Author: Ethan Coon (ecoon @ lanl.gov)

  Implementation of base class Evaluator helpers.

------------------------------------------------------------------------- */
#include "Evaluator.hh"

namespace Amanzi {
std::ostream&
operator<<(std::ostream& os, const Evaluator& self) {
  return os << self.WriteToString();
}

} // namespace
