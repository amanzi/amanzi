/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <memory>
#include <vector>

#include "errors.hh"

#include "Operator_Diagonal.hh"
#include "Op_Diagonal.hh"
#include "OperatorDefs.hh"
#include "PDE_CouplingFlux.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialize operator from parameter list.
****************************************************************** */
void PDE_CouplingFlux::Init_(
    Teuchos::ParameterList& plist,
    const Teuchos::RCP<const CompositeVectorSpace>& cvs_row,
    const Teuchos::RCP<const CompositeVectorSpace>& cvs_col,
    std::shared_ptr<const std::vector<std::vector<int> > >& row_inds,
    std::shared_ptr<const std::vector<std::vector<int> > >& col_inds)
{
  if (global_op_ == Teuchos::null) {
    global_op_ = Teuchos::rcp(new Operator_Diagonal(cvs_row, cvs_col, plist, OPERATOR_SCHEMA_INDICES));
    std::string name("Coupling_DIAGONAL");
  }

  // register the advection Op
  std::string row_compname = *(cvs_row->begin());
  std::string col_compname = *(cvs_col->begin());

  std::string name("Coupling_DIAGONAL");
  local_op_ = Teuchos::rcp(new Op_Diagonal(name, row_compname, col_compname, row_inds, col_inds));
  global_op_->OpPushBack(local_op_);
}


/* ******************************************************************
* Populate containers of elemental matrices using MFD factory.
* NOTE: input parameters are not yet used.
****************************************************************** */
void PDE_CouplingFlux::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                                      const Teuchos::Ptr<const CompositeVector>& p)
{
  auto& matrices = local_op_->matrices;
  AMANZI_ASSERT(matrices.size() == K_->size());

  WhetStone::DenseMatrix Acell;
  Acell.Reshape(1, 1);

  for (int n = 0; n < matrices.size(); ++n) {
    Acell(0, 0) = (*K_)[n] * factor_;
    matrices[n] = Acell;
  }
}

}  // namespace Operators
}  // namespace Amanzi

