/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Helper class for discrete PDE operators.
*/

#include "errors.hh"
#include "AmanziComm.hh"
#include "PDE_HelperDiscretization.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Simple constructors.
****************************************************************** */
PDE_HelperDiscretization::PDE_HelperDiscretization(const Teuchos::RCP<Operator>& global_op) :
    PDE_HelperDiscretization(global_op_->Mesh())
{
  global_op_ = global_op;
}

PDE_HelperDiscretization::PDE_HelperDiscretization(
      const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    mesh_(mesh)
{
  PopulateDimensions_();
}


/* ******************************************************************
* Supporting private routines.
****************************************************************** */
void PDE_HelperDiscretization::PopulateDimensions_()
{
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  if (mesh_->valid_edges()) {
    nedges_owned = mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::OWNED);
    nedges_wghost = mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::ALL);
  }
}


/* ******************************************************************
* Replace container of local matrices with another container.
****************************************************************** */
void PDE_HelperDiscretization::set_local_op(const Teuchos::RCP<Op>& op)
{
  if (global_op_.get()) {
    if (local_op_.get()) {
      auto index = std::find(global_op_->begin(), global_op_->end(), local_op_) - global_op_->begin();
      if (index != global_op_->size()) {
        global_op_->OpPushBack(op);
      } else {
        global_op_->OpReplace(op, index);
      }
    } else {
      global_op_->OpPushBack(op);
    }
  }
  local_op_ = op;
}

// boundary conditions (BC) require information on test and
// trial spaces. For a single PDE, these BCs could be the same.
void PDE_HelperDiscretization::SetBCs(const Teuchos::RCP<const BCs>& bc_trial,
        const Teuchos::RCP<const BCs>& bc_test)
{
  bcs_trial_.clear();
  bcs_test_.clear();
  
  bcs_trial_.emplace_back(bc_trial);
  bcs_test_.emplace_back(bc_test);
}

void PDE_HelperDiscretization::AddBCs(const Teuchos::RCP<const BCs>& bc_trial,
        const Teuchos::RCP<const BCs>& bc_test)
{
  bcs_trial_.emplace_back(bc_trial);
  bcs_test_.emplace_back(bc_test);
}


}  // namespace Operators
}  // namespace Amanzi


