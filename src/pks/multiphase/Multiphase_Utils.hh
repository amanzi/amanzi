/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Multiphase PK

*/

#ifndef AMANZI_MULTIPHASE_UTILS_HH_
#define AMANZI_MULTIPHASE_UTILS_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "Key.hh"

namespace Amanzi {
namespace Multiphase {

// no-member functions
void
ConvertFieldToTensor(const Teuchos::RCP<State>& S,
                     int dim,
                     const std::string& key,
                     std::vector<WhetStone::Tensor>& K);

KeyPair
splitPhase(const Key& name);
Key
mergePhase(const Key& name, const int phase);

} // namespace Multiphase
} // namespace Amanzi

#endif
