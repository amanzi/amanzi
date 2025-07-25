/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Collection of non-member functions that
    1. Support remapping data between different geometric objects,
       inluding coconvolution of fields such as f2 = f2 @ Map(f1).
    2. Support remapping data between different meshes.
*/

#ifndef AMANZI_OPERATORS_REMAP_UTILS_HH_
#define AMANZI_OPERATORS_REMAP_UTILS_HH_

#include "Teuchos_RCP.hpp"

namespace Amanzi {

class CompositeVector;

namespace Operators {

int CellToFace_ScaleInverse(Teuchos::RCP<const CompositeVector> f1,
                            Teuchos::RCP<CompositeVector>& f2);

} // namespace Operators
} // namespace Amanzi


#endif
