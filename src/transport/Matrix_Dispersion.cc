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

#include "Transport_constants.hh"
#include "Matrix_Dispersion.hh"


namespace Amanzi {
namespace AmanziTransport {

/* *******************************************************************
 * 
 ****************************************************************** */
void Matrix_Dispersion::Init(Dispersion_Specs& specs)
{
  specs_ = &specs;

  ncells_owned  = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  nfaces_owned  = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  dim = mesh_->space_dimension();

  // allocate memory (do it dynamically ?)
  hap_points_.resize(nfaces_wghost);
  for (int f = 0; f < nfaces_wghost; f++) hap_points_[f].init(dim);

  hap_weights_.resize(nfaces_wghost);

  D.resize(ncells_wghost);
  for (int c = 0; c < ncells_wghost; c++) D[c].init(dim, 2);
}


/* *******************************************************************
 * Calculate a dispersive tensor the from Darcy fluxes. The flux is
 * assumed to be scaled by face area.
 ****************************************************************** */
void Matrix_Dispersion::CalculateDispersionTensor(const Epetra_Vector& darcy_flux)
{
  AmanziMesh::Entity_ID_List nodes, faces;
  AmanziGeometry::Point velocity(dim), flux(dim);
  WhetStone::Tensor T(dim, 2);

  for (int c = 0; c < ncells_wghost; c++) {
    if (specs_->method == TRANSPORT_DISPERSIVITY_MODEL_ISOTROPIC) {
      for (int i = 0; i < dim; i++) D[c](i, i) = specs_->dispersivity_longitudinal;
    } else {
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      int num_good_corners = 0;
      for (int n = 0; n < nnodes; n++) {
        int v = nodes[n];
        mesh_->node_get_cell_faces(v, c, AmanziMesh::USED, &faces);
        int nfaces = faces.size();

        for (int i = 0; i < dim; i++) {
          int f = faces[i];
          const AmanziGeometry::Point& normal = mesh_->face_normal(f);
          T.add_row(i, normal);
          flux[i] = darcy_flux[f];
        }

        T.inverse();
        velocity += T * flux;
        num_good_corners++;  // each corners is good temporary (lipnikov@lanl.gov)
      }
      velocity /= num_good_corners;

      double velocity_value = norm(velocity);
      double anisotropy = specs_->dispersivity_longitudinal - specs_->dispersivity_transverse;

      for (int i = 0; i < dim; i++) {
        D[c](i, i) = specs_->dispersivity_transverse * velocity_value;
        for (int j = i; j < dim; j++) {
          double s = anisotropy * velocity[i] * velocity[j];
          if (velocity_value) s /= velocity_value;
          D[c](j, i) = D[c](i, j) += s;
        }
      }
    }

    // double vol_phi_ws = mesh_->cell_volume(c) * phi[c] * ws[c];
    // for (int i = 0; i < dim; i++) dispersion_tensor[c](i, i) *= vol_phi_ws;
  }
}


/* ******************************************************************
* Initialize Trilinos matrices. It must be called only once. 
* If matrix is non-symmetric, we generate transpose of the matrix 
* block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.   
****************************************************************** */
void Matrix_Dispersion::SymbolicAssembleGlobalMatrix()
{
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? TRANSPORT_QUAD_FACES : TRANSPORT_HEX_FACES;
  Epetra_FECrsGraph pp_graph(Copy, cmap, avg_entries_row + 1);

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
  Dpp_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, pp_graph));
  Dpp_->GlobalAssemble();
}


/* ******************************************************************
* Calculate fluxes... 
****************************************************************** */
void Matrix_Dispersion::AssembleGlobalMatrix()
{
  AmanziMesh::Entity_ID_List cells, faces;
  std::vector<int> dirs;

  // populate transmissibilities
  WhetStone::MFD3D_Diffusion mfd3d(mesh_);

  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  Epetra_Vector T(fmap_wghost);

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double> Mff(nfaces, nfaces);
    mfd3d.MassMatrixInverseTPFA(c, D[c], Mff);
   
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      T[f] += 1.0 / Mff(n, n);
    }
  }
 
  // populate the global matrix
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  Dpp_->PutScalar(0.0);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    mesh_->cell_get_face_adj_cells(c, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    // populate face-based matrix.
    Teuchos::SerialDenseMatrix<int, double> Bpp(1, ncells + 1);
    int cells_GID[ncells + 1];

    int c_GID = cmap_wghost.GID(c);  // diagonal entry
    cells_GID[0] = c_GID;

    for (int n = 0; n < ncells; n++) {
      int m = n + 1;
      cells_GID[m] = cmap_wghost.GID(cells[n]);

      int f = faces[n];
      double area = mesh_->face_area(f);
      Bpp(0, m) = -area / T[f];
      Bpp(0, 0) += area / T[f];
    }

    Dpp_->SumIntoGlobalValues(1, &c_GID, ncells + 1, cells_GID, Bpp.values());
  }
  Dpp_->GlobalAssemble();
}


/* ******************************************************************
 * * Adds time derivative to the cell-based part of MFD algebraic system.                                               
 * ****************************************************************** */
void Matrix_Dispersion::AddTimeDerivative(
    double dT, const Epetra_Vector& phi, const Epetra_Vector& ws)
{
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = volume * phi[c] * ws[c] / dT;

    int c_GID = cmap_wghost.GID(c);
    Dpp_->SumIntoGlobalValues(1, &c_GID, &factor);
  }
}


/* *******************************************************************
* Collect time-dependent boundary data in face-based arrays.                               
******************************************************************* */
void Matrix_Dispersion::Apply(const Epetra_Vector& v,  Epetra_Vector& av) const
{
  Dpp_->Apply(v, av);
}


/* *******************************************************************
* Collect time-dependent boundary data in face-based arrays.                               
******************************************************************* */
void Matrix_Dispersion::ApplyInverse(const Epetra_Vector& v,  Epetra_Vector& hv) const
{
  hv = v;
}


/* *******************************************************************
* Collect time-dependent boundary data in face-based arrays.                               
******************************************************************* */
void Matrix_Dispersion::ExtractBoundaryConditions(const int component,
                                             std::vector<int>& bc_face_id,
                                             std::vector<double>& bc_face_value)
{
  bc_face_id.assign(nfaces_wghost, 0);

  /*
  for (int n = 0; n < bcs.size(); n++) {
    if (component == bcs_tcc_index[n]) {
      for (Amanzi::Functions::TransportBoundaryFunction::Iterator bc = bcs[n]->begin(); bc != bcs[n]->end(); ++bc) {
        int f = bc->first;
        bc_face_id[f] = TRANSPORT_BC_CONSTANT_TCC;
        bc_face_value[f] = bc->second;
      }
    }
  }
  */
}


/* *******************************************************************
* Calculate harmonic averaging points and related weigths.
******************************************************************* */
void Matrix_Dispersion::PopulateHarmonicPoints()
{
  WhetStone::NLFV nlfv(mesh_);

  for (int f = 0; f < nfaces_owned; f++) {
    nlfv.HarmonicAveragingPoint(f, D, hap_points_[f], hap_weights_[f]);
  }
}


}  // namespace AmanziTransport
}  // namespace Amanzi



