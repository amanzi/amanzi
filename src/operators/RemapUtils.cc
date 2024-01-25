/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Collection of non-member functions f2 = Map(f1, f2) where 
  Map() connects fields living on different geometric objects.
*/

#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "Mesh.hh"
#include "Tensor.hh"

#include "RemapUtils.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* f2 = f2 * Map(f1)
****************************************************************** */
int CellToFace_Scale(Teuchos::RCP<CompositeVector>& f1,
                     Teuchos::RCP<CompositeVector>& f2)
{
  return 0;
}


/* ******************************************************************
* f2 = Map(f1, f2):
*   cell comp:  f2_cell = f2_cell / f1_cell
*   face comp:  f2_face = f2_face / FaceAverage(f1_cell)
****************************************************************** */
int CellToFace_ScaleInverse(Teuchos::RCP<const CompositeVector> f1,
                            Teuchos::RCP<CompositeVector>& f2)
{
  AMANZI_ASSERT(f1->hasComponent("cell"));
  AMANZI_ASSERT(f2->hasComponent("cell") && f2->hasComponent("face"));

  f1->scatterMasterToGhosted("cell");

  const Epetra_MultiVector& f1c = *f1->viewComponent("cell", true);
  Epetra_MultiVector& f2c = *f2->viewComponent("cell", true);
  Epetra_MultiVector& f2f = *f2->viewComponent("face", true);

  AmanziMesh::Entity_ID_List cells;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = f1->Map().Mesh();

  // cell-part of the map
  int ncells_wghost = mesh->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::ALL);
  for (int c = 0; c < ncells_wghost; ++c) {
    f2c[0][c] /= f1c[0][c]; 
  }

  // face-part of the map
  int nfaces_wghost = mesh->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::ALL);
  for (int f = 0; f < nfaces_wghost; ++f) {
    mesh->face_get_cells(f, AmanziMesh::Parallel_kind::ALL, &cells);
    int ncells = cells.size();

    double tmp(0.0);
    for (int n = 0; n < ncells; ++n) tmp += f1c[0][cells[n]];
    f2f[0][f] /= (tmp / ncells); 
  }

  // hack
  if (f2->hasComponent("grav")) {
    Epetra_MultiVector& f2f = *f2->viewComponent("grav", true);

    for (int f = 0; f < nfaces_wghost; ++f) {
      mesh->face_get_cells(f, AmanziMesh::Parallel_kind::ALL, &cells);
      int ncells = cells.size();

      double tmp(0.0);
      for (int n = 0; n < ncells; ++n) tmp += f1c[0][cells[n]];
      f2f[0][f] /= (tmp / ncells); 
    }
  }

  return 0;
}

}  // namespace Operators
}  // namespace Amanzi


