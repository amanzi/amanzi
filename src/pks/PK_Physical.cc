/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky
*/

/*
  Process Kernels

*/

#include "PK_Physical.hh"

#include "EvaluatorIndependentFunction.hh"
#include "State.hh"
#include "TreeVector.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Transfer operators -- copies ONLY pointers
// -----------------------------------------------------------------------------
Teuchos::RCP<TreeVectorSpace>
PK_Physical::getSolutionSpace() const
{
  CompositeVectorSpace& cvs = S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_next_);
  auto soln_space = Teuchos::rcp(new TreeVectorSpace(cvs.getComm()));
  soln_space->setData(cvs.CreateSpace());
  return soln_space;
}

void
PK_Physical::State_to_Solution(const Tag& tag, TreeVector& solution)
{
  solution.setData(S_->GetPtrW<CompositeVector>(key_, tag, getName()));
}


void
PK_Physical::Solution_to_State(const TreeVector& solution, const Tag& tag)
{
  AMANZI_ASSERT(solution.getData() == S_->GetPtr<CompositeVector>(key_, tag));
}


} // namespace Amanzi
