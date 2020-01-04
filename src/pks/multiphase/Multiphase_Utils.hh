/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MULTIPHASE_UTILS_HH_
#define AMANZI_MULTIPHASE_UTILS_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

namespace Amanzi {
namespace Multiphase {

// no-member functions
void ConvertFieldToTensor(const Teuchos::RCP<State>& S, int dim,
                          const std::string& key, std::vector<WhetStone::Tensor>& K);

}  // namespace Multiphase
}  // namespace Amanzi

#endif

