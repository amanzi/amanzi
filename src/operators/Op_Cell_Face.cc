/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#include "Op_Cell_Face.hh"

/*
  Op classes are small structs that play two roles:

  1. They provide a class name to the schema, enabling visitor patterns.
  2. They are a container for local matrices.
  
  This Op class is for storing local matrices of length ncells and with dofs
  on faces, i.e. for Schur complements in MFD methods.
*/

namespace Amanzi {
namespace Operators {

bool
Op_Cell_Face::ApplyBC(BCs& bc,
                          const Teuchos::Ptr<CompositeVector>& rhs,                       
                          bool bc_previously_applied) {

  AmanziMesh::Entity_ID_List faces;
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  const std::vector<int>& bc_model = bc.bc_model();
  const std::vector<double>& bc_value = bc.bc_value();
  const std::vector<double>& bc_mixed = bc.bc_mixed();
  ASSERT(bc_model.size() == nfaces_wghost);
  ASSERT(bc_value.size() == nfaces_wghost);

  rhs->PutScalarGhosted(0.);
  Epetra_MultiVector& rhs_face = *rhs->ViewComponent("face", true);
  Epetra_MultiVector& rhs_cell = *rhs->ViewComponent("cell");

  if (bc_previously_applied) {
    // THIS IS EXTENSIVELY UNTESTED! --etc
    ASSERT(0);
    
    for (int c=0; c!=ncells_owned; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      WhetStone::DenseMatrix& Acell = matrices[c];

      bool flag(true);
      for (int n=0; n!=nfaces; ++n) {
        int f = faces[n];
        double value = bc_value[f];

        if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
          if (flag) {  // make a copy of elemental matrix
            matrices_shadow[c] = Acell;
            flag = false;
          }
          for (int m = 0; m < nfaces; m++) {
            if (bc_model[faces[m]] != OPERATOR_BC_DIRICHLET)
              rhs_face[0][faces[m]] -= Acell(m, n) * value;
            Acell(n, m) = Acell(m, n) = 0.0;
          }
        } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
          // pass
        } else if (bc_model[f] == OPERATOR_BC_MIXED) {
          // pass
        }
      }
    }
  } else {

    for (int c=0; c!=ncells_owned; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      WhetStone::DenseMatrix& Acell = matrices[c];

      bool flag(true);
      for (int n=0; n!=nfaces; ++n) {
        int f = faces[n];
        double value = bc_value[f];

        if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
          if (flag) {  // make a copy of elemental matrix
            matrices_shadow[c] = Acell;
            flag = false;
          }
          for (int m = 0; m < nfaces; m++) {
            if (bc_model[faces[m]] != OPERATOR_BC_DIRICHLET)
              rhs_face[0][faces[m]] -= Acell(m, n) * value;
            Acell(n, m) = Acell(m, n) = 0.0;
          }
          rhs_face[0][f] = value;
          Acell(n,n) = 1.0;

        } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
          rhs_face[0][f] -= value * mesh_->face_area(f);
        } else if (bc_model[f] == OPERATOR_BC_MIXED) {
          if (flag) {  // make a copy of elemental matrix
            matrices_shadow[c] = Acell;
            flag = false;
          }
          double area = mesh_->face_area(f);
          rhs_face[0][f] -= value * area;
          Acell(n, n) += bc_mixed[f] * area;
        }
      }
    }
  }

  rhs->GatherGhostedToMaster("face", Add);
  return !bc_previously_applied;
}

}  // namespace Operators
}  // namespace Amanzi

