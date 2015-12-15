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
#include "ParallelCommunication.hh"

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
* and the positive decomposition of face conormals. The face-based
* data from left and right cells are ordered by the global cells ids.
****************************************************************** */
void OperatorDiffusionNLFV::InitStencils_()
{
  // allocate persistent memory
  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("hap", AmanziMesh::FACE, dim_)
      ->AddComponent("gamma", AmanziMesh::FACE, 1)
      ->AddComponent("weight", AmanziMesh::FACE, 2 * dim_);
  stencil_data_ = Teuchos::rcp(new CompositeVector(cvs));

  Epetra_MultiVector& hap = *stencil_data_->ViewComponent("hap", true);
  Epetra_MultiVector& gamma = *stencil_data_->ViewComponent("gamma", true);
  Epetra_MultiVector& weight = *stencil_data_->ViewComponent("weight", true);

  stencil_data_->PutScalarMasterAndGhosted(0.0);

  stencil_faces_.resize(2 * dim_);
  stencil_cells_.resize(2 * dim_);
  for (int i = 0; i < 2 * dim_; ++i) {
    stencil_faces_[i] = Teuchos::rcp(new Epetra_IntVector(mesh_->face_map(true)));
    stencil_cells_[i] = Teuchos::rcp(new Epetra_IntVector(mesh_->face_map(true)));

    stencil_faces_[i]->PutValue(0);
    stencil_cells_[i]->PutValue(0);
  }
  
  // allocate temporary memory for distributed tensor
  CompositeVectorSpace cvs_tmp; 
  cvs_tmp.SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("tensor", AmanziMesh::CELL, dim_ * dim_);
  Teuchos::RCP<CompositeVector> cv_tmp = Teuchos::rcp(new CompositeVector(cvs_tmp));
  Epetra_MultiVector& Ktmp = *cv_tmp->ViewComponent("tensor", true);

  // instantiate variables to access supporting tools
  WhetStone::NLFV nlfv(mesh_);
  WhetStone::MFD3D_Diffusion mfd3d(mesh_);

  // distribute diffusion tensor
  for (int c = 0; c < ncells_owned; ++c) {
    int k = 0;
    for (int i = 0; i < dim_; ++i) {
      for (int j = 0; j < dim_; ++j) {
        Ktmp[k][c] = (*K_)[c](i, j);
        k++;
      }
    }
  }
  cv_tmp->ScatterMasterToGhosted();

  // calculate harmonic averaging points (HAPs)
  int c1, c2;
  double hap_weight;
  AmanziMesh::Entity_ID_List cells, faces;
  AmanziGeometry::Point Kn1(dim_), Kn2(dim_), p(dim_);

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    if (ncells == 2) {
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      OrderCellsByGlobalId_(cells, c1, c2);

      Kn1 = (*K_)[c1] * normal;  // co-normals
      Kn2 = (*K_)[c2] * normal;
   
      nlfv.HarmonicAveragingPoint(f, c1, c2, Kn1, Kn2, p, hap_weight);
    } else {
      p = mesh_->face_centroid(f);
      hap_weight = 1.0;
    }

    for (int i = 0; i < dim_; ++i) hap[i][f] = p[i];
    gamma[0][f] = hap_weight;
  }

  stencil_data_->ScatterMasterToGhosted("hap");
  stencil_data_->ScatterMasterToGhosted("gamma");

  // calculate coefficients in positive decompositions of conormals
  std::vector<int> dirs;
  AmanziGeometry::Point conormal(dim_), v(dim_);
  std::vector<AmanziGeometry::Point> tau;

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    // calculate list of candidate vectors
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    tau.clear();
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      for (int i = 0; i < dim_; ++i) v[i] = hap[i][f] - xc[i];
      tau.push_back(v);
    }

    // decompose co-normals
    int ids[dim_];
    double ws[dim_];
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      conormal = ((*K_)[c] * normal) * dirs[n];

      nlfv.PositiveDecomposition(n, tau, conormal, ws, ids);

      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      OrderCellsByGlobalId_(cells, c1, c2);
      int k = (c == c1) ? 0 : dim_;

      for (int i = 0; i < dim_; i++) {
        weight[k + i][f] = ws[i];
        (*stencil_faces_[k + i])[f] = faces[ids[i]];
        (*stencil_cells_[k + i])[f] = mfd3d.cell_get_face_adj_cell(c, faces[ids[i]]);
      }
    }
  }

  // distribute stencils 
  stencil_data_->GatherGhostedToMaster("weight");
  stencil_data_->ScatterMasterToGhosted("weight");

  ParallelCommunication pp(mesh_);
  for (int i = 0; i < 2 * dim_; ++i) {
    pp.CombineGhostFace2MasterFace(*stencil_faces_[i], (Epetra_CombineMode)Add);
    pp.CombineGhostFace2MasterFace(*stencil_cells_[i], (Epetra_CombineMode)Add);

    pp.CopyMasterFace2GhostFace(*stencil_faces_[i]);
    pp.CopyMasterFace2GhostFace(*stencil_cells_[i]);
  }

  stencil_initialized_ = true;
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
  if (!stencil_initialized_) InitStencils_();

  u->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& uc = *u->ViewComponent("cell", true);

  Epetra_MultiVector& hap_gamma = *stencil_data_->ViewComponent("gamma", true);
  Epetra_MultiVector& weight = *stencil_data_->ViewComponent("weight", true);

  // allocate auxiliary matrix structure
  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_)->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, 2);
  CompositeVector matrix_cv(cvs), sideflux_cv(cvs);

  Epetra_MultiVector& matrix = *matrix_cv.ViewComponent("face", true);
  Epetra_MultiVector& sideflux = *sideflux_cv.ViewComponent("face", true);

  // calculate one-sides flux corrections. Since a flux stencil can 
  // occupy (dim_ + 1) cells, we need parallel communications.
  OneSidedFluxCorrections_(*u, sideflux_cv);

  // split each stencil between different local matrices
  int c1, c2, c3, c4, k1, k2;
  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List cells, cells_tmp, faces;

  matrix_cv.PutScalarMasterAndGhosted(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();
    
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int ncells = cells.size();

      OrderCellsByGlobalId_(cells, c1, c2);
      k1 = (c1 == c) ? 0 : 1;
      k2 = k1 * dim_;      

      // calculate solution-dependent weigths using corrections
      // to the two-point flux
      double gamma, g1, g2, gg(-1.0), w1, w2(0.0), tpfa, mu(1.0);

      gamma = hap_gamma[0][f];
      if (ncells == 2) {
        w1 = weight[0][f] * gamma;
        w2 = weight[dim_][f] * (1.0 - gamma);

        g1 = sideflux[k1][f];
        g2 = sideflux[1 - k1][f];
        gg = g1 * g2;

        g1 = fabs(g1);
        g2 = fabs(g2);
        mu = (g1 + g2 == 0.0) ? 0.5 : g2 / (g1 + g2);
      } else {
        w1 = weight[0][f];
      }

      tpfa = mu * w1 + (1.0 - mu) * w2;
      matrix[k1][f] += tpfa;

      // remaining terms of one-sided flux in cell c
      if (gg <= 0.0) {
        for (int i = 1; i < dim_; i++) {
          int f1 = (*stencil_faces_[i + k2])[f];
          mesh_->face_get_cells(f1, AmanziMesh::USED, &cells_tmp);

          gamma = hap_gamma[0][f1];
          OrderCellsByGlobalId_(cells_tmp, c3, c4);

          k1 = 0;
          if (c3 != c) {
            gamma = 1.0 - gamma;
            k1 = 1;
          }

          double tmp = ncells * weight[i + k2][f] * gamma * mu;
          matrix[k1][f1] += tmp;
        }
      }
    }
  }

  matrix_cv.GatherGhostedToMaster();

  // populate local matrices
  for (int f = 0; f < nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    WhetStone::DenseMatrix Aface(ncells, ncells);

    if (ncells == 2) {
      int k1 = OrderCellsByGlobalId_(cells, c3, c4);
      int k2 = 1 - k1;
      Aface(0, 0) = matrix[k1][f];
      Aface(0, 1) = -matrix[k1][f];

      Aface(1, 0) = -matrix[k2][f];
      Aface(1, 1) = matrix[k2][f];
    } else {
      Aface(0, 0) = matrix[0][f];
    }
  
    local_op_->matrices[f] = Aface;
  }
}


/* ******************************************************************
* Calculate one-sided flux for given cell c and face f.
****************************************************************** */
double OperatorDiffusionNLFV::OneSidedFluxCorrections_(
  const CompositeVector& u, CompositeVector& flux_cv) 
{
  const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);
  Epetra_MultiVector& flux = *flux_cv.ViewComponent("face", true);

  Epetra_MultiVector& hap_gamma = *stencil_data_->ViewComponent("gamma", true);
  Epetra_MultiVector& weight = *stencil_data_->ViewComponent("weight", true);

  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();
  
  int c1, c2, c3, k1, k2;
  double gamma, tmp;
  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List cells, cells_tmp, faces;

  flux_cv.PutScalarMasterAndGhosted(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();
    
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);

      OrderCellsByGlobalId_(cells, c1, c2);
      k1 = (c1 == c) ? 0 : 1;
      k2 = k1 * dim_;      
    
      double sideflux(0.0);
      for (int i = 1; i < dim_; ++i) {
        int f1 = (*stencil_faces_[i + k2])[f];
        c3 = (*stencil_cells_[i + k2])[f];
        if (c3 >= 0) {
          mesh_->face_get_cells(f1, AmanziMesh::USED, &cells_tmp);
          OrderCellsByGlobalId_(cells_tmp, c1, c2);

          gamma = hap_gamma[0][f1];
          if (c1 != c) gamma = 1.0 - gamma;

          tmp = weight[i + k2][f] * gamma;
          sideflux += tmp * (uc[0][c] - uc[0][c3]);
        } else {
          tmp = weight[i + k2][f];
          sideflux += tmp * (uc[0][c] - bc_value[f1]);
        }
      }
      flux[k1][f] = sideflux; 
    }
  }

  flux_cv.GatherGhostedToMaster();
  flux_cv.ScatterMasterToGhosted();
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


/* ******************************************************************
* Order cells by their global ids. Returns 1 if cells were swapped.
****************************************************************** */
int OperatorDiffusionNLFV::OrderCellsByGlobalId_(
    const AmanziMesh::Entity_ID_List& cells, int& c1, int& c2)
{
  c1 = cells[0];
  c2 = -1;

  int ncells = cells.size();
  if (ncells == 1) return 0;

  c2 = cells[1];
  if (mesh_->cell_map(true).GID(c1) > mesh_->cell_map(true).GID(c2)) {
    int c(c1);
    c1 = c2;
    c2 = c;
    return 1;
  } 

  return 0;
}

}  // namespace Operators
}  // namespace Amanzi

