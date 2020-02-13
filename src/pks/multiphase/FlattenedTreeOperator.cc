/*
  Multiphase PK

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Schema.hh"

#include "FlattenedTreeOperator.hh"

namespace Amanzi {
namespace Multiphase {

FlattenedTreeOperator::FlattenedTreeOperator(Teuchos::RCP<const TreeVectorSpace> tvs) 
{
  block_diagonal_ = false;

  Operators::Schema schema;
  tvs_flat_ = Teuchos::rcp(new TreeVectorSpace());

  int n_blocks(0);
  for (auto it = tvs->begin(); it != tvs->end(); ++it) {
    auto cvs = (*it)->Data();
    int nvec = cvs->NumVectors(*cvs->begin());

    auto cvs0 = Teuchos::rcp(new CompositeVectorSpace());
    cvs0->SetMesh(cvs->Mesh())->SetGhosted(true);
    cvs0->SetOwned(false);

    for (auto kt = cvs->begin(); kt != cvs->end(); ++kt) {
      AMANZI_ASSERT(nvec == cvs->NumVectors(*kt));
      cvs0->AddComponent(*kt, schema.StringToKind(*kt), 1);  // FIXME
    }
    CompositeVector cv(*cvs0);

    for (int n = 0; n < nvec; ++n) {
      Teuchos::RCP<TreeVectorSpace> tvs0 = Teuchos::rcp(new TreeVectorSpace());
      tvs0->SetData(cvs0);
      tvs_flat_->PushBack(tvs0);
    }
    n_blocks += nvec;
  }

  // resize the blocks
  blocks_.resize(n_blocks, Teuchos::Array<Teuchos::RCP<const Operators::Operator> >(n_blocks, Teuchos::null));
  tvs_ = tvs_flat_;
  // tvs_->Print(std::cout);
}


/* ******************************************************************
* Calculate Y = A * X using matrix-free matvec on blocks of operators.
****************************************************************** */
int FlattenedTreeOperator::Apply(const TreeVector& X, TreeVector& Y) const {
  return ApplyAssembled(X, Y);
}

}  // namespace Multiphase
}  // namespace Amanzi
