/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MECHANICS_FRACTURED_MATRIX_PK_HH_
#define AMANZI_MECHANICS_FRACTURED_MATRIX_PK_HH_

#include <string>
#include <vector>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

// Amanzi::Mechanics
#include "MechanicsSmallStrain_PK.hh"

namespace Amanzi {
namespace Mechanics {

class MechanicsFracturedMatrix_PK : public MechanicsSmallStrain_PK {
 public:
  MechanicsFracturedMatrix_PK(Teuchos::ParameterList& pk_tree,
                              const Teuchos::RCP<Teuchos::ParameterList>& glist,
                              const Teuchos::RCP<State>& S,
                              const Teuchos::RCP<TreeVector>& soln);

  ~MechanicsFracturedMatrix_PK(){};
};

} // namespace Mechanics
} // namespace Amanzi

#endif

