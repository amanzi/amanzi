/*
  State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Implementation of base class Evaluator helpers.
*/

#include "Evaluator.hh"

namespace Amanzi {
std::ostream &operator<<(std::ostream &os, const Evaluator& self) {
  return os << self.WriteToString();
}

}  // namespace Amanzi
