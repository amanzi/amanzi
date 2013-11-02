/*
This is the transport component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Teuchos_RCP.hpp"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "mfd3d_diffusion.hh"
#include "nlfv.hh"
#include "tensor.hh"
#include "PreconditionerFactory.hh"

#include "TransportDefs.hh"
#include "Dispersion_TPFA.hh"

namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
* Initialize Trilinos matrices. It must be called only once. 
* If matrix is non-symmetric, we generate transpose of the matrix 
* block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.   
****************************************************************** */
void Dispersion_TPFA::SymbolicAssembleMatrix()
{
  const Epetra_Map& cmap_owned = mesh_->cell_map(false);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  int avg_entries_row = (dim == 2) ? TRANSPORT_QUAD_FACES : TRANSPORT_HEX_FACES;
  Epetra_FECrsGraph pp_graph(Copy, cmap_owned, avg_entries_row + 1);

  AmanziMesh::Entity_ID_List cells;
  int cells_GID[2];

  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    for (int n = 0; n < ncells; n++)
        cells_GID[n] = cmap_wghost.GID(cells[n]);

    pp_graph.InsertGlobalIndices(ncells, cells_GID, ncells, cells_GID);
  }
  pp_graph.GlobalAssemble();  // Symbolic graph is complete.

  // create global matrices
  App_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, pp_graph));
  App_->GlobalAssemble();
}


/* ******************************************************************
* Calculate and assemble fluxes using the TPFA scheme.
****************************************************************** */
void Dispersion_TPFA::AssembleMatrix(const Epetra_MultiVector& p)
{
  AmanziMesh::Entity_ID_List cells, faces;
  std::vector<int> dirs;

  // populate transmissibilities
  WhetStone::MFD3D_Diffusion mfd3d(mesh_);

  Teuchos::RCP<CompositeVector> T = CreateCompositeVector(mesh_, AmanziMesh::FACE, 1, true);
  T->CreateData();
  Teuchos::RCP<Epetra_MultiVector> Ttmp = T->ViewComponent("face", true);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    WhetStone::DenseMatrix Mff(nfaces, nfaces);
    mfd3d.MassMatrixInverseTPFA(c, D[c], Mff);
   
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      (*Ttmp)[0][f] += 1.0 / Mff(n, n);
    }
  }
  T->GatherGhostedToMaster();
 
  // populate the global matrix
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  int cells_GID[2];
  WhetStone::DenseMatrix Bpp(2, 2);

  App_->PutScalar(0.0);

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    if (ncells < 2) continue;

    for (int n = 0; n < ncells; n++) {
      cells_GID[n] = cmap_wghost.GID(cells[n]);

      double coef = 1.0 / (*Ttmp)[0][f];
      Bpp(0, 0) =  coef;
      Bpp(1, 1) =  coef;
      Bpp(0, 1) = -coef;
      Bpp(1, 0) = -coef;
    }

    App_->SumIntoGlobalValues(ncells, cells_GID, Bpp.Values());
  }
  App_->GlobalAssemble();
}


/* *******************************************************************
* Collect time-dependent boundary data in face-based arrays.                               
******************************************************************* */
void Dispersion_TPFA::Apply(const Epetra_Vector& v, Epetra_Vector& av) const
{
  App_->Apply(v, av);
}


/* *******************************************************************
* Collect time-dependent boundary data in face-based arrays.                               
******************************************************************* */
int Dispersion_TPFA::ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const
{
  preconditioner_->ApplyInverse(v, hv);
  return 0;
}


}  // namespace AmanziTransport
}  // namespace Amanzi



