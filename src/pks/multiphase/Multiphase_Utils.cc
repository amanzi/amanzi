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

#include <string>
#include <vector>

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "State.hh"
#include "Tensor.hh"

#include "MultiphaseDefs.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Convertion from double to tensor.
****************************************************************** */
void
ConvertFieldToTensor(const Teuchos::RCP<State>& S,
                     int dim,
                     const std::string& key,
                     std::vector<WhetStone::Tensor>& K)
{
  const auto& cv = S->Get<CompositeVector>(key);
  cv.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& perm = *cv.ViewComponent("cell", true);

  int ncells = perm.MyLength();
  int ndofs = perm.NumVectors();
  K.resize(ncells);
  bool off_diag = cv.HasComponent("offd");

  // most common cases of diagonal permeability
  if (ndofs == 1) {
    for (int c = 0; c < ncells; c++) {
      K[c].Init(dim, 1);
      K[c](0, 0) = perm[0][c];
    }
  } else if (dim == 2) {
    for (int c = 0; c < ncells; c++) {
      if (!off_diag && perm[0][c] == perm[1][c]) {
        K[c].Init(dim, 1);
        K[c](0, 0) = perm[0][c];
      } else {
        K[c].Init(dim, 2);
        K[c](0, 0) = perm[0][c];
        K[c](1, 1) = perm[1][c];
      }
    }
  } else if (dim == 3) {
    for (int c = 0; c < K.size(); c++) {
      if (!off_diag && perm[0][c] == perm[1][c] && perm[0][c] == perm[2][c]) {
        K[c].Init(dim, 1);
        K[c](0, 0) = perm[0][c];
      } else {
        K[c].Init(dim, 2);
        K[c](0, 0) = perm[0][c];
        K[c](1, 1) = perm[1][c];
        K[c](2, 2) = perm[2][c];
      }
    }
  }

  // special case of permeability with off-diagonal components
  if (off_diag) {
    const Epetra_MultiVector& offd = *cv.ViewComponent("offd");

    for (int c = 0; c < ncells; c++) {
      if (dim == 2) {
        K[c](0, 1) = K[c](1, 0) = offd[0][c];
      } else if (dim == 3) {
        K[c](0, 1) = K[c](1, 0) = offd[0][c];
        K[c](0, 2) = K[c](2, 0) = offd[1][c];
        K[c](1, 2) = K[c](2, 1) = offd[2][c];
      }
    }
  }
}


/* ******************************************************************
* Parsing phase suffix
****************************************************************** */
KeyPair
splitPhase(const Key& name)
{
  std::size_t pos = name.find_last_of(Keys::phase_delimiter);
  if (pos == std::string::npos) return std::make_pair(Key(""), name);
  else return std::make_pair(name.substr(0, pos), name.substr(pos + 1, std::string::npos));
}


Key
mergePhase(const Key& name, const int phase)
{
  if (phase == MULTIPHASE_PHASE_LIQUID) return name + Keys::phase_delimiter + "liquid";
  else if (phase == MULTIPHASE_PHASE_GAS) return name + Keys::phase_delimiter + "gas";
  else return name + Keys::phase_delimiter + "napl";
}

} // namespace Multiphase
} // namespace Amanzi
