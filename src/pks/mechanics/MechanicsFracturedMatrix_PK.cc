/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Mechanics PK

*/

#include <vector>

#include "MechanicsFracturedMatrix_PK.hh"

namespace Amanzi {
namespace Mechanics {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* New constructor: extracts lists and requires fields.
****************************************************************** */
MechanicsFracturedMatrix_PK::MechanicsFracturedMatrix_PK(Teuchos::ParameterList& pk_tree,
                                                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                                         const Teuchos::RCP<State>& S,
                                                         const Teuchos::RCP<TreeVector>& soln)
  : MechanicsSmallStrain_PK(pk_tree, glist, S, soln)
{
}

} // namespace Mechanics
} // namespace Amanzi

