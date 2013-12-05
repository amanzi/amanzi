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
#include "DenseMatrix.hh"
#include "PreconditionerFactory.hh"

#include "TransportDefs.hh"
#include "Dispersion_MFD.hh"

namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
* Initialize Trilinos matrices. It must be called only once. 
* If matrix is non-symmetric, we generate transpose of the matrix 
* block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.   
****************************************************************** */
void Dispersion_MFD::SymbolicAssembleMatrix()
{
  // create supermap
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);

  int ncells_global = cmap.NumGlobalElements();
  int nfaces_global = fmap.NumGlobalElements();

  int ndof_local = ncells_owned + nfaces_owned;
  int ndof_global = ncells_global + nfaces_global; 

  int* gids = new int[ndof_local];
  cmap.MyGlobalElements(&(gids[0]));
  fmap.MyGlobalElements(&(gids[ncells_owned]));

  for (int f = 0; f < nfaces_owned; f++) gids[ncells_owned + f] += ncells_global;
  super_map_ = Teuchos::rcp(new Epetra_Map(ndof_global, ndof_local, gids, 0, cmap.Comm()));

  delete [] gids;

  // create graph of the matrix sparcity structure
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  int avg_entries_row = (dim == 2) ? TRANSPORT_QUAD_FACES : TRANSPORT_HEX_FACES;
  Epetra_FECrsGraph graph(Copy, *super_map_, 2 * avg_entries_row);

  AmanziMesh::Entity_ID_List faces, cells;
  std::vector<int> dirs;
  int dof_GID[TRANSPORT_MAX_FACES + 1];

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();
    int ndof = nfaces + 1;

    dof_GID[0] = cmap.GID(c);
    for (int n = 0; n < nfaces; n++) {
      dof_GID[n + 1] = ncells_global + fmap_wghost.GID(faces[n]);
    }
  
    graph.InsertGlobalIndices(ndof, dof_GID, ndof, dof_GID);
  }
  graph.GlobalAssemble();

  // create matrix
  App_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, graph));
  App_->GlobalAssemble();
}


/* ******************************************************************
* Calculate and assemble fluxes using the MFD method.
****************************************************************** */
void Dispersion_MFD::AssembleMatrix(const Epetra_Vector& p)
{
  WhetStone::MFD3D_Diffusion mfd3d(mesh_);

  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  int ncells_global = cmap.NumGlobalElements();
  int dof_GID[TRANSPORT_MAX_FACES + 1];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  App_->PutScalar(0.0);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();
    int ndof = nfaces + 1;

    dof_GID[0] = cmap.GID(c);
    for (int n = 0; n < nfaces; n++) {
      dof_GID[n + 1] = ncells_global + fmap_wghost.GID(faces[n]);
    }

    WhetStone::DenseMatrix Bpp(ndof, ndof);
    WhetStone::DenseMatrix Wff(nfaces, nfaces);

    // build inverse of the mass matrix
    int i = mfd3d.MassMatrixInverseMMatrixHex(c, D[c], Wff);
    if (i == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_WRONG) {
      mfd3d.MassMatrixInverseOptimizedScaled(c, D[c], Wff);
    }
    for (int n = 0; n < nfaces; n++)
      for (int m = 0; m < nfaces; m++) Bpp(m + 1, n + 1) = Wff(m, n);

    double matsum = 0.0;
    for (int n = 0; n < nfaces; n++) {
      double rowsum = 0.0;
      for (int m = 0; m < nfaces; m++) rowsum += Wff(n, m);

      Bpp(n + 1, 0) = Bpp(0, n + 1) = -rowsum;
      matsum += rowsum;
    }
    Bpp(0, 0) = matsum;

    App_->SumIntoGlobalValues(ndof, dof_GID, Bpp.Values());
  }

  App_->GlobalAssemble();
}


/* *******************************************************************
* Collect time-dependent boundary data in face-based arrays.                               
******************************************************************* */
int Dispersion_MFD::Apply(const Epetra_Vector& v, Epetra_Vector& av) const
{
  return App_->Apply(v, av);
}


/* *******************************************************************
* Collect time-dependent boundary data in face-based arrays.                               
******************************************************************* */
int Dispersion_MFD::ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const
{
  preconditioner_->ApplyInverse(v, hv);
  return 0;
}


}  // namespace AmanziTransport
}  // namespace Amanzi



