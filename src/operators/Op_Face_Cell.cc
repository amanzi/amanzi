/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#include "Op_Face_Cell.hh"

/*
  Op classes are small structs that play two roles:

  1. They provide a class name to the schema, enabling visitor patterns.
  2. They are a container for local matrices.
  
  This Op class is for storing local matrices of length nfaces and with dofs
  on cells, i.e. for Advection or for TPFA.
*/

namespace Amanzi {
namespace Operators {

bool
Op_Face_Cell::ApplyBC(BCs& bc,
                      const Teuchos::Ptr<CompositeVector>& rhs,                       
                      bool bc_previously_applied) {

  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  AmanziMesh::Entity_ID_List cells, nodes;

  const std::vector<int>& bc_model = bc.bc_model();
  const std::vector<double>& bc_value = bc.bc_value();
  const std::vector<double>& bc_mixed = bc.bc_mixed();
  ASSERT(bc_model.size() == nfaces_wghost);
  ASSERT(bc_value.size() == nfaces_wghost);

  Epetra_MultiVector& rhs_cell = *rhs->ViewComponent("cell");

  for (int f = 0; f != nfaces_owned; ++f) {
    WhetStone::DenseMatrix& Aface = matrices[f];

    if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      rhs_cell[0][cells[0]] += bc_value[f] * Aface(0, 0);
    }
    else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
      matrices_shadow[f] = Aface;

      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      rhs_cell[0][cells[0]] -= bc_value[f] * mesh_->face_area(f);
      Aface *= 0.0;
    }
    // solve system of two equations in three unknowns
    else if (bc_model[f] == OPERATOR_BC_MIXED) {
      matrices_shadow[f] = Aface;

      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      double area = mesh_->face_area(f);
      double factor = area / (1.0 + bc_mixed[f] * area / Aface(0, 0));
      rhs_cell[0][cells[0]] -= bc_value[f] * factor;
      Aface(0, 0) = bc_mixed[f] * factor;
    }
  }

  return false; // we do not directly apply a face-based BC
}


}  // namespace Operators
}  // namespace Amanzi


