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
void Matrix_Dispersion::CalculateDispersionTensor(const Epetra_Vector& darcy_flux, 
                                                  const Epetra_Vector& porosity, 
                                                  const Epetra_Vector& saturation)
{
  if (specs_->model == TRANSPORT_DISPERSIVITY_MODEL_ISOTROPIC) {
    for (int c = 0; c < ncells_wghost; c++) {
      for (int i = 0; i < dim; i++) D[c](i, i) = specs_->dispersivity_longitudinal;
    }
  } else {
    WhetStone::MFD3D_Diffusion mfd3d(mesh_);

    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;
    AmanziGeometry::Point velocity(dim);

    for (int c = 0; c < ncells_wghost; c++) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      std::vector<double> flux(nfaces);
      for (int n = 0; n < nfaces; n++) flux[n] = darcy_flux[faces[n]];
      mfd3d.RecoverGradient_MassMatrix(c, flux, velocity);

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

      D[c] *= porosity[c] * saturation[c];
    }
  }
}


/* ******************************************************************
* Initialize Trilinos matrices. It must be called only once. 
* If matrix is non-symmetric, we generate transpose of the matrix 
* block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.   
****************************************************************** */
void Matrix_Dispersion::SymbolicAssembleGlobalMatrix()
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
  int cells_GID[2];
  Teuchos::SerialDenseMatrix<int, double> Bpp(2, 2);

  Dpp_->PutScalar(0.0);

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    if (ncells < 2) continue;

    for (int n = 0; n < ncells; n++) {
      cells_GID[n] = cmap_wghost.GID(cells[n]);

      double coef = mesh_->face_area(f) / T[f];
      Bpp(0, 0) =  coef;
      Bpp(1, 1) =  coef;
      Bpp(0, 1) = -coef;
      Bpp(1, 0) = -coef;
    }

    Dpp_->SumIntoGlobalValues(ncells, cells_GID, Bpp.values());
  }
  Dpp_->GlobalAssemble();
}


/* ******************************************************************
* Adds time derivative to the cell-based part of MFD algebraic system.
****************************************************************** */
void Matrix_Dispersion::AddTimeDerivative(
    double dT, const Epetra_Vector& porosity, const Epetra_Vector& saturation)
{
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = volume * porosity[c] * saturation[c] / dT;

    int c_GID = cmap_wghost.GID(c);
    Dpp_->SumIntoGlobalValues(1, &c_GID, &factor);
  }
}


/* *******************************************************************
* Collect time-dependent boundary data in face-based arrays.                               
******************************************************************* */
void Matrix_Dispersion::Apply(const Epetra_Vector& v, Epetra_Vector& av) const
{
  Dpp_->Apply(v, av);
}


/* *******************************************************************
* Collect time-dependent boundary data in face-based arrays.                               
******************************************************************* */
void Matrix_Dispersion::ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const
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



