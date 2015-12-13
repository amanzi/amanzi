/*
   Operators

   Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
   Amanzi is released under the three-clause BSD License. 
   The terms of use and "as is" disclaimer for this license are 
   provided in the top-level COPYRIGHT file.

   Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
            Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// Amanzi
#include "mfd3d_diffusion.hh"
#include "nlfv.hh"

#include "OperatorDefs.hh"
#include "OperatorDiffusionNLFV.hh"
#include "Op_Face_Cell.hh"
#include "Operator_Cell.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization
****************************************************************** */
void OperatorDiffusionNLFV::InitDiffusion_(Teuchos::ParameterList& plist)
{
  // Define stencil for the FV diffusion method.
  local_op_schema_ = OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL;

  // create or check the existing Operator
  if (global_op_ == Teuchos::null) {
    // constructor was given a mesh
    global_op_schema_ = OPERATOR_SCHEMA_DOFS_CELL;

    // build the CVS from the global schema
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);
    cvs->AddComponent("cell", AmanziMesh::CELL, 1);

    global_op_ = Teuchos::rcp(new Operator_Cell(cvs, plist, global_op_schema_));

  } else {
    // constructor was given an Operator
    global_op_schema_ = global_op_->schema();
    mesh_ = global_op_->DomainMap().Mesh();
  }

  // create the local Op and register it with the global Operator
  std::string name = "Diffusion: FACE_CELL";
  local_op_ = Teuchos::rcp(new Op_Face_Cell(name, mesh_));
  global_op_->OpPushBack(local_op_);

  // mesh info
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  // other data
  dim_ = mesh_->space_dimension();
}


/* ******************************************************************
* Compute harmonic averaging points (function of geometry and tensor)
* and the positive decomposition of face conormals.
****************************************************************** */
void OperatorDiffusionNLFV::InitStencils_()
{
  // allocate memory
  stencils_.resize(nfaces_owned);
  for (int f = 0; f < nfaces_owned; f++) stencils_[f].Init(dim_);
  
  WhetStone::NLFV nlfv(mesh_);
  WhetStone::MFD3D_Diffusion mfd3d(mesh_);

  // calculate harmonic averaging points (HAPs)
  for (int f = 0; f < nfaces_owned; f++) {
    nlfv.HarmonicAveragingPoint(f, *K_, stencils_[f].p, stencils_[f].gamma);
  }

  // calculate coefficients in positive decompositions of conormals
  AmanziMesh::Entity_ID_List cells, faces;
  std::vector<int> dirs;

  AmanziGeometry::Point conormal(dim_), v(dim_);
  std::vector<AmanziGeometry::Point> tau;

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    // calculate local directions from centroid to HAPs
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    tau.clear();
    for (int i = 0; i < nfaces; i++) {
      v = stencils_[faces[i]].p - xc;
      tau.push_back(v);
    }

    // calculate positive decomposition of exterior conormals
    int ids[dim_];
    double ws[dim_];
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      if (f < nfaces_owned) {
        const AmanziGeometry::Point& normal = mesh_->face_normal(f);
        conormal = ((*K_)[c] * normal) * dirs[n];

        nlfv.PositiveDecomposition(n, tau, conormal, ws, ids);

        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        int k = (cells[0] == c) ? 0 : dim_;

        for (int i = 0; i < dim_; i++) {
          stencils_[f].weights[k + i] = ws[i];
          stencils_[f].stencil[k + i] = mfd3d.cell_get_face_adj_cell(c, faces[ids[i]]);
          stencils_[f].faces[k + i] = faces[ids[i]];
        }
      }
    }
  }

  stencils_initialized_ = true;
}


/* ******************************************************************
* Populate face-based 2x2 matrices on interior faces and 1x1 matrices
* on boundary faces. We avoid round-off operations since the stencils 
* already incorporate them.
****************************************************************** */
void OperatorDiffusionNLFV::UpdateMatrices(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u)
{
  if (!stencils_initialized_) InitStencils_();

  u->ScatterMasterToGhosted("cell");

  // allocate zero local matrices
  AmanziMesh::Entity_ID_List cells, cells_tmp;

  for (int f = 0; f < nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    WhetStone::DenseMatrix Aface(ncells, ncells);
    Aface = 0.0;
    local_op_->matrices[f] = Aface;
  }
 
  // split each stencil between different local matrices
  const Epetra_MultiVector& uc = *u->ViewComponent("cell", true);

  for (int f = 0; f < nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    // calculate solution-dependent weigths using corrections
    // to the two-point flux
    int c1, c2, f1;
    double gamma, tmp, g1, g2, gg(-1.0), w1, w2(0.0), tpfa, mu(1.0);

    gamma = stencils_[f].gamma;
    c1 = cells[0];

    if (ncells == 2) {
      c2 = cells[1];
      w1 = stencils_[f].weights[0] * gamma;
      w2 = stencils_[f].weights[dim_] * (1.0 - gamma);

      g1 = OneSidedFluxCorrection_(c1, f, uc, 0);
      g2 = OneSidedFluxCorrection_(c2, f, uc, dim_);
      gg = g1 * g2;

      g1 = fabs(g1);
      g2 = fabs(g2);
      mu = (g1 + g2 == 0.0) ? 0.5 : g2 / (g1 + g2);
    } else {
      w1 = stencils_[f].weights[0];
    }
    tpfa = mu * w1 + (1.0 - mu) * w2;

    // add the TPFA term of the flux to the local matrix
    WhetStone::DenseMatrix& Aface = local_op_->matrices[f];
    for (int i = 0; i < ncells; ++i) {
      Aface(i, i) += tpfa;
      for (int j = i + 1; j < ncells; ++j) {
        Aface(i, j) -= tpfa;
        Aface(j, i) -= tpfa;
      }
    }

    if (gg <= 0.0) {
      // calculate the remaining terms of the left flux
      for (int i = 1; i < dim_; i++) {
        f1 = stencils_[f].faces[i];
        if (f1 >= nfaces_owned) continue;
        mesh_->face_get_cells(f1, AmanziMesh::USED, &cells_tmp);

        int k(0);
        gamma = stencils_[f1].gamma;
        if (cells_tmp[0] != c1) {
          gamma = 1.0 - gamma;
          k = 1;
        }

        tmp = ncells * stencils_[f].weights[i] * gamma * mu;
        WhetStone::DenseMatrix& Bface = local_op_->matrices[f1];
        Bface(k, k) += tmp;
        if (Bface.NumRows() == 2) Bface(k, 1 - k) -= tmp;
      }
    
      // calculate the remaining terms of the right flux
      for (int i = 1; i < dim_; i++) {
        f1 = stencils_[f].faces[i + dim_];
        if (f1 >= nfaces_owned) continue;
        mesh_->face_get_cells(f1, AmanziMesh::USED, &cells_tmp);

        int k(0);
        gamma = stencils_[f1].gamma;
        if (cells_tmp[0] != c2) {
          gamma = 1.0 - gamma;
          k = 1;
        } 

        tmp = ncells * stencils_[f].weights[i + dim_] * gamma * (1.0 - mu);
        WhetStone::DenseMatrix& Bface = local_op_->matrices[f1];
        Bface(k, k) += tmp;
        if (Bface.NumRows() == 2) Bface(k, 1 - k) -= tmp;
      }
    }
  }
}


/* ******************************************************************
* Calculate one-sided flux for given cell c and face f.
****************************************************************** */
double OperatorDiffusionNLFV::OneSidedFluxCorrection_(
    int c, int f, const Epetra_MultiVector& uc, int k) 
{
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();
  
  int c3, f1;
  double gamma, tmp, flux(0.0);
  AmanziMesh::Entity_ID_List cells;

  for (int i = 1; i < dim_; i++) {
    f1 = stencils_[f].faces[i + k];
    if (f1 >= nfaces_owned) continue;
    mesh_->face_get_cells(f1, AmanziMesh::USED, &cells);

    c3 = stencils_[f].stencil[i + k];
    if (c3 >= 0) {
      gamma = stencils_[f1].gamma;
      if (cells[0] != c) gamma = 1.0 - gamma;

      tmp = stencils_[f].weights[i + k] * gamma;
      flux += tmp * (uc[0][c] - uc[0][c3]);
    } else {
      tmp = stencils_[f].weights[i + k];
      flux += tmp * (uc[0][c] - bc_value[f1]); 
    }
  }

  return flux;
}


/* ******************************************************************
* Matrix-based implementation of boundary conditions.
****************************************************************** */
void OperatorDiffusionNLFV::ApplyBCs(bool primary, bool eliminate)
{
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell", true);

  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] != OPERATOR_BC_NONE) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int c = cells[0];

      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        WhetStone::DenseMatrix& Aface = local_op_->matrices[f];
        rhs_cell[0][c] += Aface(0, 0) * bc_value[f];
      } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
        local_op_->matrices_shadow[f] = local_op_->matrices[f];
        local_op_->matrices[f](0,0) = 0.0;
            
        rhs_cell[0][c] -= bc_value[f] * mesh_->face_area(f);
      }
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi

