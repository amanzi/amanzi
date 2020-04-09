/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATORS_REMAP_UTILS_HH_
#define AMANZI_OPERATORS_REMAP_UTILS_HH_

#include "Teuchos_RCP.hpp"

namespace Amanzi {

class CompositeVector;

namespace AmanziMesh {
class Mesh;
}

namespace Operators {

int
CellToFace_Scale(Teuchos::RCP<CompositeVector> f1,
                 Teuchos::RCP<CompositeVector>& f2);
int
CellToFace_ScaleInverse(Teuchos::RCP<const CompositeVector> f1,
                        Teuchos::RCP<CompositeVector>& f2);

} // namespace Operators
} // namespace Amanzi


#endif
