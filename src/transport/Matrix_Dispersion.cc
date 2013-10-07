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
#include "Matrix_Dispersion.hh"

namespace Amanzi {
namespace AmanziTransport {

/* *******************************************************************
* Initialization of a class.
******************************************************************* */
void Matrix_Dispersion::Init(std::vector<Teuchos::RCP<DispersionModel> >& specs, 
                             const std::string& preconditioner, 
                             const Teuchos::ParameterList& prec_list)
{
  specs_ = &specs;

  ncells_owned  = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  nfaces_owned  = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  dim = mesh_->space_dimension();

  D.resize(ncells_owned);
  for (int c = 0; c < ncells_owned; c++) D[c].init(dim, 2);

  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(preconditioner, prec_list);
}


/* *******************************************************************
 * Calculate a dispersive tensor the from Darcy fluxes. The flux is
 * assumed to be scaled by face area.
 ****************************************************************** */
void Matrix_Dispersion::CalculateDispersionTensor(const Epetra_Vector& darcy_flux, 
                                                  const Epetra_Vector& porosity, 
                                                  const Epetra_Vector& saturation)
{
  for (int mb = 0; mb < specs_->size(); mb++) {
    Teuchos::RCP<DispersionModel> spec = (*specs_)[mb]; 

    std::vector<AmanziMesh::Entity_ID> block;
    for (int r = 0; r < (spec->regions).size(); r++) {
      std::string region = (spec->regions)[r];
      mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

      AmanziMesh::Entity_ID_List::iterator c;
      for (c = block.begin(); c != block.end(); c++) {
        D[*c].PutScalar(0.0); 
        if (spec->model == TRANSPORT_DISPERSIVITY_MODEL_ISOTROPIC) {
          for (int i = 0; i < dim; i++) {
            D[*c](i, i) = spec->alphaL + spec->D * spec->tau;
          }
          D[*c] *= porosity[*c] * saturation[*c];
        } else {
          WhetStone::MFD3D_Diffusion mfd3d(mesh_);

          AmanziMesh::Entity_ID_List faces;
          std::vector<int> dirs;
          AmanziGeometry::Point velocity(dim);

          mesh_->cell_get_faces_and_dirs(*c, &faces, &dirs);
          int nfaces = faces.size();

          std::vector<double> flux(nfaces);
          for (int n = 0; n < nfaces; n++) flux[n] = darcy_flux[faces[n]];
          mfd3d.RecoverGradient_MassMatrix(*c, flux, velocity);
          velocity /= porosity[*c];  // pore velocity

          double velocity_value = norm(velocity);
          double anisotropy = spec->alphaL - spec->alphaT;

          for (int i = 0; i < dim; i++) {
            D[*c](i, i) = spec->D * spec->tau + spec->alphaT * velocity_value;
            for (int j = i; j < dim; j++) {
              double s = anisotropy * velocity[i] * velocity[j];
              if (velocity_value) s /= velocity_value;
              D[*c](j, i) = D[*c](i, j) += s;
            }
          }

          D[*c] *= porosity[*c] * saturation[*c];
        }
      }
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
  App_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, pp_graph));
  App_->GlobalAssemble();
}


/* ******************************************************************
* Calculate and assemble fluxes using the TPFA scheme.
****************************************************************** */
void Matrix_Dispersion::AssembleGlobalMatrixTPFA(const Teuchos::RCP<Transport_State>& TS)
{
  AmanziMesh::Entity_ID_List cells, faces;
  std::vector<int> dirs;

  // populate transmissibilities
  WhetStone::MFD3D_Diffusion mfd3d(mesh_);

  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  Epetra_Vector T(fmap_wghost);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    WhetStone::DenseMatrix Mff(nfaces, nfaces);
    mfd3d.MassMatrixInverseTPFA(c, D[c], Mff);
   
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      T[f] += 1.0 / Mff(n, n);
    }
  }
  TS->CombineGhostFace2MasterFace(T, Add);
 
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

      double coef = 1.0 / T[f];
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
* Allocate necessary structures for the nonlinear scheme.
******************************************************************* */
void Matrix_Dispersion::InitNLFV()
{
  int d = mesh_->space_dimension();
  stencil_.resize(nfaces_wghost);
  for (int f = 0; f < nfaces_wghost; f++) stencil_[f].Init(d);
}


/* ******************************************************************
* Create face-based flux stencils.
****************************************************************** */
void Matrix_Dispersion::CreateFluxStencils()
{
  WhetStone::NLFV nlfv(mesh_);
  WhetStone::MFD3D_Diffusion mfd3d(mesh_);

  // calculate harmonic averaging points
  for (int f = 0; f < nfaces_owned; f++) {
    nlfv.HarmonicAveragingPoint(f, D, stencil_[f].p, stencil_[f].gamma);
  }

  // calculate coefficients in positive decompositions of conormals
  AmanziMesh::Entity_ID_List cells, faces;
  std::vector<int> dirs;

  int d = mesh_->space_dimension();
  AmanziGeometry::Point conormal(d), v(d);
  std::vector<AmanziGeometry::Point> tau;

  std::vector<int> fpointer;
  fpointer.assign(nfaces_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    // calculate local directions
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    tau.clear();
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      v = stencil_[f].p - xc;
      tau.push_back(v);
    }

    // calculate positive decomposition of the conormals
    int ids[d];
    double ws[d];
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      conormal = (D[c] * normal) * dirs[n];

      nlfv.PositiveDecomposition(n, tau, conormal, ws, ids);

      int k = fpointer[f];
      for (int i = 0; i < d; i++) {
        stencil_[f].weights[k + i] = ws[i];
        stencil_[f].stencil[k + i] = mfd3d.cell_get_face_adj_cell(c, faces[ids[i]]);
        stencil_[f].faces[k + i] = faces[ids[i]];
      }
      fpointer[f] += d;
    }
  }
}


/* ******************************************************************
* Calculate and assemble fluxes using the NLFV scheme. We avoid 
* round-off operations since the stencils already incorporate them.
****************************************************************** */
void Matrix_Dispersion::AssembleGlobalMatrixNLFV(const Epetra_Vector& p)
{
  AmanziMesh::Entity_ID_List cells;
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

  int d = mesh_->space_dimension();
  int cells_GID[d + 1];
  double Bpp[d + 1];

  // populate the global matrix (loop over internal faces)
  App_->PutScalar(0.0);

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    if (ncells < 2) continue;

    // Select cell with the smallest id. Since \gamma is 
    // associated with the first cell it must be swapped.
    double gamma_tpfa = stencil_[f].gamma;
    int c1 = cells[0];
    int c2 = cells[1];
    if (c1 > c2) {
      int c = c1;
      c1 = c2;
      c2 = c;
      gamma_tpfa = 1.0 - gamma_tpfa;
    }

    // calculate non-TPFA fluxes. Typically g1 * g2 <= 0.0
    double g1 = 0.0;
    for (int i = 1; i < d; i++) {
      int c3 = stencil_[f].stencil[i];
      if (c3 >= 0) {
        int f1 = stencil_[f].faces[i];
        double gamma = stencil_[f1].gamma;
        if (c3 < c1) gamma = 1.0 - gamma;

        double tmp = stencil_[f].weights[i] *gamma;
        g1 += tmp * (p[c1] - p[c3]);
      }
    }

    double g2 = 0.0;
    for (int i = 1; i < d; i++) {
      int c3 = stencil_[f].stencil[i + d];
      if (c3 >= 0) {
        int f1 = stencil_[f].faces[i + d];
        double gamma = stencil_[f1].gamma;
        if (c3 < c2) gamma = 1.0 - gamma;

        double tmp = stencil_[f].weights[i + d] * gamma;
        g2 += tmp * (p[c2] - p[c3]);
      }
    }

    // calculate TPFA flux
    double w1, w2, tpfa, gg, mu(0.5);
    gg = g1 * g2;
    g1 = fabs(g1);
    g2 = fabs(g2);
    if (g1 + g2 != 0.0) mu = g2 / (g1 + g2);

    w1 = stencil_[f].weights[0] * gamma_tpfa;
    w2 = stencil_[f].weights[d] * gamma_tpfa;
    tpfa = mu * w1 + (1.0 - mu) * w2;

    // add fluxes to the matrix
    c1 = cmap_wghost.GID(c1);
    c2 = cmap_wghost.GID(c2);

    int m = 1;
    cells_GID[0] = c1;
    cells_GID[1] = c2;
    Bpp[0] =  tpfa;
    Bpp[1] = -tpfa;

    if (gg <= 0.0) { 
      for (int i = 1; i < d; i++) {
        int c3 = stencil_[f].stencil[i];
        if (c3 >= 0) {
          m++;
          cells_GID[m] = cmap_wghost.GID(c3);

          int f1 = stencil_[f].faces[i];
          double gamma = stencil_[f1].gamma;
          if (c3 < c1) gamma = 1.0 - gamma;

          double tmp = 2 * stencil_[f].weights[i] * gamma * mu;
          Bpp[0] += tmp;
          Bpp[m] = -tmp;
        }
      }
    }

    App_->SumIntoGlobalValues(1, &c1, m + 1, cells_GID, Bpp);

    m = 1;
    cells_GID[0] = c2;
    cells_GID[1] = c1;
    Bpp[0] =  tpfa;
    Bpp[1] = -tpfa;

    if (gg <= 0.0) { 
      for (int i = 1; i < d; i++) {
        int c3 = stencil_[f].stencil[i + d];
        if (c3 >= 0) {
          m++;
          cells_GID[m] = cmap_wghost.GID(c3);

          int f1 = stencil_[f].faces[i + d];
          double gamma = stencil_[f1].gamma;
          if (c3 < c2) gamma = 1.0 - gamma;

          double tmp = 2 * stencil_[f].weights[i + d] * gamma * (1.0 - mu);
          Bpp[0] += tmp;
          Bpp[m] = -tmp;
        }
      }
    }

    App_->SumIntoGlobalValues(1, &c2, m + 1, cells_GID, Bpp);
  }
  App_->GlobalAssemble();
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
    App_->SumIntoGlobalValues(1, &c_GID, &factor);
  }
}


/* *******************************************************************
* Collect time-dependent boundary data in face-based arrays.                               
******************************************************************* */
void Matrix_Dispersion::Apply(const Epetra_Vector& v, Epetra_Vector& av) const
{
  App_->Apply(v, av);
}


/* *******************************************************************
* Collect time-dependent boundary data in face-based arrays.                               
******************************************************************* */
void Matrix_Dispersion::ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const
{
  preconditioner_->ApplyInverse(v, hv);
}


}  // namespace AmanziTransport
}  // namespace Amanzi



