/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "Epetra_FECrsGraph.h"

#include "mfd3d_diffusion.hh"
#include "PreconditionerFactory.hh"

#include "Flow_PK.hh"
#include "FlowDefs.hh"
#include "Matrix_MFD.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Constructor                                      
****************************************************************** */
Matrix_MFD::Matrix_MFD(Teuchos::RCP<State> S,
                       std::vector<WhetStone::Tensor>* K, 
                       Teuchos::RCP<RelativePermeability> rel_perm)
    : Matrix<CompositeVector, CompositeVectorSpace>(S, K, rel_perm)
{ 
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  actions_ = 0;

  npassed_ = 0;
  nokay_ = 0;
}


/* ******************************************************************
* Cannot destroy ML cleanly. Try Trilinos 10.10 (lipnikov@lanl.gov)
****************************************************************** */
Matrix_MFD::~Matrix_MFD()
{
  // if (MLprec->IsPreconditionerComputed()) {
  //   MLprec->DestroyPreconditioner();
  //   delete MLprec;
  // }
}


/* ******************************************************************
* Calculate elemental inverse mass matrices. 
* WARNING: The original Aff matrices are destroyed.                                            
****************************************************************** */
void Matrix_MFD::CreateMassMatrices(int mfd3d_method)
{
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);

  Mff_cells_.clear();

  int ok;
  nokay_ = npassed_ = 0;

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    int nfaces = mesh_->cell_get_num_faces(c);

    WhetStone::DenseMatrix Mff(nfaces, nfaces);
    WhetStone::Tensor& Kc = (*K_)[c];

    if (mfd3d_method == FLOW_MFD3D_POLYHEDRA_SCALED) {
      ok = mfd.MassMatrixInverseScaled(c, Kc, Mff);
    } else if (mfd3d_method == FLOW_MFD3D_POLYHEDRA_MONOTONE) {
      ok = mfd.MassMatrixInverseMMatrix(c, Kc, Mff);
      if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_WRONG) {
        ok = mfd.MassMatrixInverseTPFA(c, Kc, Mff);
        nokay_--;
        npassed_++;
      } 
    } else if (mfd3d_method == FLOW_MFD3D_POLYHEDRA) {
      ok = mfd.MassMatrixInverse(c, Kc, Mff);
    } else if (mfd3d_method == FLOW_MFD3D_OPTIMIZED_SCALED) {
      ok = mfd.MassMatrixInverseOptimizedScaled(c, Kc, Mff);
    } else if (mfd3d_method == FLOW_MFD3D_OPTIMIZED) {
      ok = mfd.MassMatrixInverseOptimized(c, Kc, Mff);
    } else if (mfd3d_method == FLOW_MFD3D_HEXAHEDRA_MONOTONE) {
      if ((nfaces == 6 && dim == 3) || (nfaces == 4 && dim == 2))
        ok = mfd.MassMatrixInverseMMatrixHex(c, Kc, Mff);
      else
        ok = mfd.MassMatrixInverse(c, Kc, Mff);
    } else if (mfd3d_method == FLOW_MFD3D_TPFA) {
      ok = mfd.MassMatrixInverseTPFA(c, Kc, Mff);
    } else if (mfd3d_method == FLOW_MFD3D_SUPPORT_OPERATOR) {
      ok = mfd.MassMatrixInverseSO(c, Kc, Mff);
    } else {
      Errors::Message msg("Flow PK: unexpected discretization methods (contact lipnikov@lanl.gov).");
      Exceptions::amanzi_throw(msg);
    }

    Mff_cells_.push_back(Mff);

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
      Errors::Message msg("Flow PK: unexpected failure of LAPACK in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }
    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_OK) nokay_++;
    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_PASSED) npassed_++;
  }

  // sum up the numbers across processors
  int nokay_tmp = nokay_, npassed_tmp = npassed_;
  mesh_->get_comm()->SumAll(&nokay_tmp, &nokay_, 1);
  mesh_->get_comm()->SumAll(&npassed_tmp, &npassed_, 1);
}


/* ******************************************************************
* Calculate elemental inverse mass matrices.                                           
****************************************************************** */
void Matrix_MFD::CreateMassMatrices_ScaledStability(int mfd3d_method, double factor)
{
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);

  mfd.ModifyStabilityScalingFactor(factor);

  Mff_cells_.clear();

  int ok;
  nokay_ = npassed_ = 0;

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    int nfaces = mesh_->cell_get_num_faces(c);

    WhetStone::DenseMatrix Mff(nfaces, nfaces);
    WhetStone::Tensor& Kc = (*K_)[c];

    if (mfd3d_method == FLOW_MFD3D_HEXAHEDRA_MONOTONE) {
      if ((nfaces == 6 && dim == 3) || (nfaces == 4 && dim == 2))
        ok = mfd.MassMatrixInverseMMatrixHex(c, Kc, Mff);
      else
        ok = mfd.MassMatrixInverse(c, Kc, Mff);
    } else if (mfd3d_method == FLOW_MFD3D_POLYHEDRA_MONOTONE) {
      ok = mfd.MassMatrixInverseMMatrix(c, Kc, Mff);
      if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_WRONG) {
        ok = mfd.MassMatrixInverseTPFA(c, Kc, Mff);
      }
    } else if (mfd3d_method == FLOW_MFD3D_TPFA) {
      ok = mfd.MassMatrixInverseTPFA(c, Kc, Mff);
    } else if (mfd3d_method == FLOW_MFD3D_SUPPORT_OPERATOR) {
      ok = mfd.MassMatrixInverseSO(c, Kc, Mff);
    } else if (mfd3d_method == FLOW_MFD3D_OPTIMIZED) {
      ok = mfd.MassMatrixInverseOptimized(c, Kc, Mff);
    } else {
      ok = mfd.MassMatrixInverse(c, Kc, Mff);
    }

    Mff_cells_.push_back(Mff);

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
      Errors::Message msg("Matrix_MFD: unexpected failure of LAPACK in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }
    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_OK) nokay_++;
    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_PASSED) npassed_++;
  }

  // sum up the numbers across processors
  int nokay_tmp = nokay_, npassed_tmp = npassed_;
  mesh_->get_comm()->SumAll(&nokay_tmp, &nokay_, 1);
  mesh_->get_comm()->SumAll(&npassed_tmp, &npassed_, 1);
}


/* ******************************************************************
* Calculate elemental stiffness matrices (fully saturated flow)                                          
****************************************************************** */
void Matrix_MFD::CreateStiffnessMatricesDarcy(int mfd3d_method)
{
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;

  Aff_cells_.clear();
  Afc_cells_.clear();
  Acf_cells_.clear();
  Acc_cells_.clear();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    int nfaces = mesh_->cell_get_num_faces(c);

    WhetStone::DenseMatrix& Mff = Mff_cells_[c];
    Teuchos::SerialDenseMatrix<int, double> Bff(nfaces, nfaces);
    Epetra_SerialDenseVector Bcf(nfaces);

    double* braw = Bff.values();
    double* mraw = Mff.Values();
    for (int n = 0; n < nfaces * nfaces; n++) braw[n] = mraw[n];

    double matsum = 0.0;  // elimination of mass matrix
    for (int n = 0; n < nfaces; n++) {
      double rowsum = 0.0;
      for (int m = 0; m < nfaces; m++) rowsum += Bff(n, m);

      Bcf(n) = -rowsum;
      matsum += rowsum;
    }

    Aff_cells_.push_back(Bff);  // This the only place where memory can be allocated.
    Afc_cells_.push_back(Bcf);
    Acf_cells_.push_back(Bcf);
    Acc_cells_.push_back(matsum);
  }
}


/* ******************************************************************
* Calculate stiffness matrices for partially saturated flow.
* Permeability class carries information on nonlinear coefficient.
****************************************************************** */
void Matrix_MFD::CreateStiffnessMatricesRichards()
{
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;

  Aff_cells_.clear();
  Afc_cells_.clear();
  Acf_cells_.clear();
  Acc_cells_.clear();

  Epetra_MultiVector& Krel_cells = *rel_perm_->Krel().ViewComponent("cell");
  Epetra_MultiVector& Krel_faces = *rel_perm_->Krel().ViewComponent("face", true);
  int method = rel_perm_->method();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    int nfaces = mesh_->cell_get_num_faces(c);

    WhetStone::DenseMatrix& Mff = Mff_cells_[c];
    Teuchos::SerialDenseMatrix<int, double> Bff(nfaces, nfaces);
    Epetra_SerialDenseVector Bcf(nfaces), Bfc(nfaces);

    if (method == FLOW_RELATIVE_PERM_NONE) {
      double* braw = Bff.values();
      double* mraw = Mff.Values();
      for (int n = 0; n < nfaces * nfaces; n++) braw[n] = mraw[n];

    } else if (method == FLOW_RELATIVE_PERM_CENTERED) {  // centered permeability for diffusion
      for (int n = 0; n < nfaces; n++)
        for (int m = 0; m < nfaces; m++) Bff(m, n) = Mff(m, n) * Krel_cells[0][c];

    } else if (method == FLOW_RELATIVE_PERM_AMANZI) {
      std::vector<double>& krel = rel_perm_->Krel_amanzi()[c];

      double tmp = Krel_cells[0][c];
      for (int n = 0; n < nfaces; n++)
        for (int m = 0; m < nfaces; m++) Bff(m, n) = Mff(m, n) * tmp;

      // add upwind correction
      mesh_->cell_get_faces(c,&faces);
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        // double t = std::max(0.0, Krel_faces[0][f] - tmp);
        double t = fabs(Krel_faces[0][f] - tmp);
        Bff(n, n) += Mff(n, n) * t; 
      }
    } else {
      mesh_->cell_get_faces(c,&faces);
      for (int m = 0; m < nfaces; m++) {
        AmanziMesh::Entity_ID f = faces[m];
        double tmp = Krel_faces[0][f];
        for (int n = 0; n < nfaces; n++) Bff(m, n) = Mff(m, n) * tmp;
      }
    }

    double matsum = 0.0;  // elimination of mass matrix
    for (int n = 0; n < nfaces; n++) {
      double rowsum = 0.0, colsum = 0.0;
      for (int m = 0; m < nfaces; m++) {
        colsum += Bff(m, n);
        rowsum += Bff(n, m);
      }
      Bcf(n) = -colsum;
      Bfc(n) = -rowsum;
      matsum += colsum;
    }

    Aff_cells_.push_back(Bff);  // This the only place where memory can be allocated.
    Afc_cells_.push_back(Bfc);
    Acf_cells_.push_back(Bcf);
    Acc_cells_.push_back(matsum);
  }
}


/* ******************************************************************
* May be used in the future.                                            
****************************************************************** */
void Matrix_MFD::RescaleMFDstiffnessMatrices(const Epetra_Vector& old_scale,
                                             const Epetra_Vector& new_scale)
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells; c++) {
    Teuchos::SerialDenseMatrix<int, double>& Bff = Aff_cells_[c];
    Epetra_SerialDenseVector& Bcf = Acf_cells_[c];

    int n = Bff.numRows();
    double scale = old_scale[c] / new_scale[c];

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) Bff(i, j) *= scale;
      Bcf(i) *= scale;
    }
    Acc_cells_[c] *= scale;
  }
}


/* ******************************************************************
* Simply allocates memory.                                           
****************************************************************** */
void Matrix_MFD::CreateRHSVectors()
{
  Ff_cells_.clear();
  Fc_cells_.clear();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;

  for (int c = 0; c < ncells; c++) {
    int nfaces = mesh_->cell_get_num_faces(c);

    Epetra_SerialDenseVector Ff(nfaces);  // Entries are initilaized to 0.0.
    double Fc = 0.0;

    Ff_cells_.push_back(Ff);
    Fc_cells_.push_back(Fc);
  }
}


/* ******************************************************************
* Routine updates elemental discretization matrices and must be 
* called before applying boundary conditions and global assembling. 
* simplified implementation for single phase flow.                                            
****************************************************************** */
void Matrix_MFD::AddGravityFluxesDarcy(double rho, const AmanziGeometry::Point& gravity)
{
  AmanziGeometry::Point rho_gravity(gravity);
  rho_gravity *= rho;

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Epetra_SerialDenseVector& Ff = Ff_cells_[c];
    double& Fc = Fc_cells_[c];

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);

      double outward_flux = (((*K_)[c] * rho_gravity) * normal) * dirs[n]; 
      Ff[n] += outward_flux;
      Fc -= outward_flux;
    }
  }
}


/* ******************************************************************
* Routine updates elemental discretization matrices and must be 
* called before applying boundary conditions and global assembling.                                             
****************************************************************** */
void Matrix_MFD::AddGravityFluxesRichards(double rho, const AmanziGeometry::Point& gravity,
                                          std::vector<int>& bc_model) 
{
  AmanziGeometry::Point rho_gravity(gravity);
  rho_gravity *= rho;

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Epetra_MultiVector& Krel_cells = *rel_perm_->Krel().ViewComponent("cell");
  Epetra_MultiVector& Krel_faces = *rel_perm_->Krel().ViewComponent("face", true);
  int method = rel_perm_->method();

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Epetra_SerialDenseVector& Ff = Ff_cells_[c];
    double& Fc = Fc_cells_[c];
    std::vector<double>& krel = rel_perm_->Krel_amanzi()[c];

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);

      double outward_flux = (((*K_)[c] * rho_gravity) * normal) * dirs[n]; 
      if (method == FLOW_RELATIVE_PERM_CENTERED) {
        outward_flux *= Krel_cells[0][c];
      } else if (method == FLOW_RELATIVE_PERM_AMANZI) {
        outward_flux *= krel[n]; 
      } else {
        outward_flux *= Krel_faces[0][f];
      }
      Ff[n] += outward_flux;
      Fc -= outward_flux;  // Nonzero-sum contribution when flag_upwind = false.
    }
  }
}


/* ******************************************************************
* Applies boundary conditions to elemental stiffness matrices and
* creates elemental rigth-hand-sides.                                           
****************************************************************** */
void Matrix_MFD::ApplyBoundaryConditions(
    std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values)
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double>& Bff = Aff_cells_[c];  // B means elemental.
    Epetra_SerialDenseVector& Bfc = Afc_cells_[c];
    Epetra_SerialDenseVector& Bcf = Acf_cells_[c];

    Epetra_SerialDenseVector& Ff = Ff_cells_[c];
    double& Fc = Fc_cells_[c];

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      double value = bc_values[f][0];

      if (bc_model[f] == FLOW_BC_FACE_PRESSURE) {
        for (int m = 0; m < nfaces; m++) {
          Ff[m] -= Bff(m, n) * value;
          Bff(n, m) = Bff(m, n) = 0.0;
        }
        Fc -= Bcf(n) * value;
        Bcf(n) = Bfc(n) = 0.0;

        Bff(n, n) = 1.0;
        Ff[n] = value;
      } else if (bc_model[f] == FLOW_BC_FACE_FLUX) {
        Ff[n] -= value * mesh_->face_area(f);
      } else if (bc_model[f] == FLOW_BC_FACE_MIXED) {
        double area = mesh_->face_area(f);
        Ff[n] += value * area;
        Bff(n, n) -= bc_values[f][1] * area;
      }

      // If one wants to deposit infiltration in soil.
      // if (bc_model[f] == FLOW_BC_FACE_PRESSURE_SEEPAGE) {
      //   Fc -= bc_values[f][1] * mesh_->face_area(f);
      // }
    }
  }
}


/* ******************************************************************
* Adds time derivative to the cell-based part of MFD algebraic system.
****************************************************************** */
void Matrix_MFD::AddTimeDerivative(
    const Epetra_MultiVector& p, const Epetra_MultiVector& phi, double rho, double dT)
{
  Epetra_MultiVector dSdP(p);
  rel_perm_->DerivedSdP(p, dSdP);

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = rho * phi[0][c] * dSdP[0][c] * volume / dT;
    Acc_cells_[c] += factor;
    Fc_cells_[c] += factor * p[0][c];
  }
}


/* ******************************************************************
* Adds time derivative related to specific stroage to cell-based 
* part of MFD algebraic system. 
****************************************************************** */
void Matrix_MFD::AddTimeDerivativeSpecificStorage(
    const Epetra_MultiVector& p, const Epetra_MultiVector& ss, double g, double dT)
{
  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = volume * ss[0][c] / (g * dT);
    Acc_cells_[c] += factor;
    Fc_cells_[c] += factor * p[0][c];
  }
}



/* ******************************************************************
* Adds time derivative related to specific yiled to cell-based part 
* of MFD algebraic system. Area factor is alreafy inside Sy. 
****************************************************************** */
void Matrix_MFD::AddTimeDerivativeSpecificYield(
    const Epetra_MultiVector& p, const Epetra_MultiVector& sy, double g, double dT)
{
  for (int c = 0; c < ncells_owned; c++) {
    double factor = sy[0][c] / (g * dT);
    Acc_cells_[c] += factor;
    Fc_cells_[c] += factor * p[0][c];
  }
}


/* ******************************************************************
* Initialize Trilinos matrices. It must be called only once. 
* If matrix is non-symmetric, we generate transpose of the matrix 
* block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.   
****************************************************************** */
void Matrix_MFD::SymbolicAssemble()
{
  // create the p-lambda map associated with the matrix
  if (cvs_.size() == 0) {  // ugly solution (lipnikov@lanl.gov) 
    cvs_.SetMesh(mesh_);
    cvs_.SetGhosted(true);
    cvs_.SetComponent("cell", AmanziMesh::CELL, 1);
    cvs_.SetOwned(false);
    cvs_.AddComponent("face", AmanziMesh::FACE, 1);
  }

  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? FLOW_QUAD_FACES : FLOW_HEX_FACES;
  Epetra_CrsGraph cf_graph(Copy, cmap, fmap_wghost, avg_entries_row, false);  // FIX (lipnikov@lanl.gov)
  Epetra_FECrsGraph ff_graph(Copy, fmap, 2*avg_entries_row);

  AmanziMesh::Entity_ID_List faces;
  int faces_GID[FLOW_MAX_FACES];

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      faces_GID[n] = fmap_wghost.GID(faces[n]);
    }
    cf_graph.InsertMyIndices(c, nfaces, &(faces[0]));
    ff_graph.InsertGlobalIndices(nfaces, faces_GID, nfaces, faces_GID);
  }
  cf_graph.FillComplete(fmap, cmap);
  ff_graph.GlobalAssemble();  // Symbolic graph is complete.

  // create global matrices
  Acc_ = Teuchos::rcp(new Epetra_Vector(cmap));
  Acf_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cf_graph));
  if (actions_ & AmanziFlow::FLOW_MATRIX_ACTION_MATRIX) {
    Aff_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
    Aff_->GlobalAssemble();
  } 
  if (actions_ & AmanziFlow::FLOW_MATRIX_ACTION_PRECONDITIONER) {
    Sff_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
    Sff_->GlobalAssemble();
  }

  if (flag_symmetry_) {
    Afc_ = Acf_;
  } else {
    Afc_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cf_graph));
  }

  rhs_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));
}


/* ******************************************************************
* Assemble elemental mass matrices into four global matrices. 
* We need an auxiliary GHOST-based vector to assemble the RHS.
****************************************************************** */
void Matrix_MFD::Assemble()
{
  Aff_->PutScalar(0.0);

  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  AmanziMesh::Entity_ID_List faces;
  int faces_GID[FLOW_MAX_FACES];

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  Epetra_MultiVector& rhs_cells = *rhs_->ViewComponent("cell");
  Epetra_MultiVector& rhs_faces = *rhs_->ViewComponent("face", true);

  rhs_faces.PutScalar(0.0);

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      faces_GID[n] = fmap_wghost.GID(faces[n]);
    }
    (*Acc_)[c] = Acc_cells_[c];
    Acf_->ReplaceMyValues(c, nfaces, Acf_cells_[c].Values(), &(faces[0]));
    Aff_->SumIntoGlobalValues(nfaces, faces_GID, Aff_cells_[c].values());

    if (!flag_symmetry_)
      Afc_->ReplaceMyValues(c, nfaces, Afc_cells_[c].Values(), &(faces[0]));

    rhs_cells[0][c] = Fc_cells_[c];
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      rhs_faces[0][f] += Ff_cells_[c][n];
    }
  }
  Aff_->GlobalAssemble();

  rhs_->GatherGhostedToMaster("face", Add);
}


/* ******************************************************************
* Assembles four matrices: diagonal Acc_, two off-diagonal blocks
* Acf_ and Afc_, and the Schur complement Sff_.
****************************************************************** */
void Matrix_MFD::AssembleSchurComplement_(
    std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values)
{
  Sff_->PutScalar(0.0);

  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  AmanziMesh::Entity_ID_List faces;
  int faces_GID[FLOW_MAX_FACES];

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();
    Teuchos::SerialDenseMatrix<int, double> Schur(nfaces, nfaces);

    Epetra_SerialDenseVector& Bcf = Acf_cells_[c];
    Epetra_SerialDenseVector& Bfc = Afc_cells_[c];

    for (int n = 0; n < nfaces; n++) {
      for (int m = 0; m < nfaces; m++) {
        Schur(n, m) = Aff_cells_[c](n, m) - Bfc[n] * Bcf[m] / Acc_cells_[c];
      }
    }

    for (int n = 0; n < nfaces; n++) {  // Symbolic boundary conditions
      int f = faces[n];
      if (bc_model[f] == FLOW_BC_FACE_PRESSURE) {
        for (int m = 0; m < nfaces; m++) Schur(n, m) = Schur(m, n) = 0.0;
        Schur(n, n) = 1.0;
      }
    }

    for (int n = 0; n < nfaces; n++) {
      faces_GID[n] = fmap_wghost.GID(faces[n]);
    }
    (*Sff_).SumIntoGlobalValues(nfaces, faces_GID, Schur.values());

    // check that the other matrices were not calculated already
    if (!(actions_ & AmanziFlow::FLOW_MATRIX_ACTION_MATRIX)) {
      (*Acc_)[c] = Acc_cells_[c];
      Acf_->ReplaceMyValues(c, nfaces, Acf_cells_[c].Values(), &(faces[0]));

      if (!flag_symmetry_)
        Afc_->ReplaceMyValues(c, nfaces, Afc_cells_[c].Values(), &(faces[0]));
    }
  }
  (*Sff_).GlobalAssemble();
}


/* ******************************************************************
* Initialization of the preconditioner                                                 
****************************************************************** */
void Matrix_MFD::InitPreconditioner(const std::string& prec_name, const Teuchos::ParameterList& prec_list)
{
  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(prec_name, prec_list);
}


/* ******************************************************************
* Parallel matvec product A * X.                                              
****************************************************************** */
int Matrix_MFD::Apply(const CompositeVector& X, CompositeVector& Y) const
{
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");
  const Epetra_MultiVector& Xf = *X.ViewComponent("face");

  Epetra_MultiVector& Yc = *Y.ViewComponent("cell");
  Epetra_MultiVector& Yf = *Y.ViewComponent("face");

  // Face unknowns:  Yf = Aff * Xf + Afc * Xc
  int ierr;
  Epetra_MultiVector Tf(Xf.Map(), 1);
  ierr  = (*Aff_).Multiply(false, Xf, Yf);
  ierr |= (*Afc_).Multiply(true, Xc, Tf);  // Afc is kept in the transpose form
  Yf.Update(1.0, Tf, 1.0);

  // Cell unknowns:  Yc = Acf * Xf + Acc * Xc
  ierr |= (*Acf_).Multiply(false, Xf, Yc);  // It performs the required parallel communications.
  ierr |= Yc.Multiply(1.0, *Acc_, Xc, 1.0);

  if (ierr) {
    Errors::Message msg("Matrix_MFD::Apply has failed to calculate y = A*x.");
    Exceptions::amanzi_throw(msg);
  }
  return 0;
}


/* ******************************************************************
* The OWNED cell-based and face-based d.o.f. are packed together into 
* the X and Y Epetra vectors, with the cell-based in the first part.
****************************************************************** */
int Matrix_MFD::ApplyPreconditioner(const CompositeVector& X, CompositeVector& Y) const
{
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");
  const Epetra_MultiVector& Xf = *X.ViewComponent("face");

  Epetra_MultiVector& Yc = *Y.ViewComponent("cell");
  Epetra_MultiVector& Yf = *Y.ViewComponent("face");

  // Temporary cell and face vectors.
  Epetra_MultiVector Tc(Xc);
  Epetra_MultiVector Tf(Xf);

  // FORWARD ELIMINATION:  Tf = Xf - Afc inv(Acc) Xc
  int ierr;
  ierr  = Tc.ReciprocalMultiply(1.0, *Acc_, Xc, 0.0);
  ierr |= (*Afc_).Multiply(true, Tc, Tf);  // Afc is kept in transpose form
  Tf.Update(1.0, Xf, -1.0);

  // Solve the Schur complement system Sff * Yf = Tf.
  preconditioner_->ApplyInverse(Tf, Yf);

  // BACKWARD SUBSTITUTION:  Yc = inv(Acc) (Xc - Acf Yf)
  ierr |= (*Acf_).Multiply(false, Yf, Tc);  // It performs the required parallel communications.
  Tc.Update(1.0, Xc, -1.0);
  ierr |= Yc.ReciprocalMultiply(1.0, *Acc_, Tc, 0.0);

  if (ierr) {
    Errors::Message msg("Matrix_MFD::ApplyInverse has failed in calculating y = inv(A)*x.");
    Exceptions::amanzi_throw(msg);
  }
  return 0;
}


/* ******************************************************************
* Reduce the pressure-lambda-system to lambda-system via ellimination
* of the known pressure. Structure of the global system is preserved
* but off-diagoal blocks are zeroed-out.                                               
****************************************************************** */
int Matrix_MFD::ReduceGlobalSystem2LambdaSystem(CompositeVector& u)
{
  Epetra_MultiVector& uc = *u.ViewComponent("cell");
  Epetra_MultiVector& gf = *rhs_->ViewComponent("face");

  // Update RHS: rhs = rhs - Afc * uc
  Epetra_MultiVector tf(gf.Map(), 1);
  Afc_->Multiply(true, uc, tf);  // Afc is kept in the transpose form.
  gf.Update(-1.0, tf, 1.0);

  // Decouple pressure-lambda system
  Afc_->PutScalar(0.0);
  Acf_->PutScalar(0.0);

  return 0;
}


/* ******************************************************************
* Replaces block of preconditioner by blocks of matrix.                                               
****************************************************************** */
int Matrix_MFD::PopulatePreconditioner(Matrix_MFD& matrix)
{
  if (actions_ & AmanziFlow::FLOW_MATRIX_ACTION_PRECONDITIONER && 
      matrix.CheckActionProperty(AmanziFlow::FLOW_MATRIX_ACTION_MATRIX)) {
    *Sff_ = *(matrix.Aff());
    *Acf_ = *(matrix.Acf());
    *Afc_ = *(matrix.Afc());
    *Acc_ = *(matrix.Acc());
    return 1;
  } else {
    return 0;
  }
}


/* ******************************************************************
* WARNING: Routines requires original mass matrices (Aff_cells), i.e.
* before boundary conditions were imposed.
*
* WARNING: Since diffusive flux is not continuous, we derive it only
* once (using flag) and in exactly the same manner as in routine
* Flow_PK::addGravityFluxes_DarcyFlux.
****************************************************************** */
void Matrix_MFD::DeriveMassFlux(
    const CompositeVector& solution, CompositeVector& darcy_mass_flux,
    std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values)
{
  solution.ScatterMasterToGhosted("face");

  const Epetra_MultiVector& solution_cells = *solution.ViewComponent("cell");
  const Epetra_MultiVector& solution_faces = *solution.ViewComponent("face", true);
  Epetra_MultiVector& flux = *darcy_mass_flux.ViewComponent("face", true);

  AmanziMesh::Entity_ID_List faces;
  std::vector<double> dp;
  std::vector<int> dirs;
  std::vector<int> flag(nfaces_wghost, 0);


  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    dp.resize(nfaces);
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      dp[n] = solution_cells[0][c] - solution_faces[0][f];
    }

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
  
      if (f < nfaces_owned && !flag[f]) {
        double s(0.0);
// if (f==11){
//   std::cout<<"number "<< n<<"\n";
//   for (int m = 0; m < nfaces; m++){	    
//     std::cout<<Aff_cells_[c](n, m) <<" "<< dp[m]<<" cell "<<solution_cells[0][c]<<" face "<<solution_faces[0][faces[m]]<<std::endl;
//   }
// }
        for (int m = 0; m < nfaces; m++) s += Aff_cells_[c](n, m) * dp[m];
        flux[0][f] = s * dirs[n];  
//if (f==11) std::cout<<"FFFF "<<flux[0][f]<<std::endl;
        flag[f] = 1;
      }
    }
  }

  if (&*rel_perm_){
    AddGravityFluxes_DarcyFlux(flux, *rel_perm_);
  }
  else {
    AddGravityFluxes_DarcyFlux(flux);
  }

}


/* ******************************************************************
* Updates global Darcy vector calculated by a discretization method.
* simplified implementation for a single phase flow.                                            
****************************************************************** */
void Matrix_MFD::AddGravityFluxes_DarcyFlux(Epetra_MultiVector& darcy_mass_flux)
{
  double rho = *(S_->GetScalarData("fluid_density"));
  const Epetra_Vector& gravity_ = *S_->GetConstantVectorData("gravity");
  int dim = mesh_->space_dimension();

  AmanziGeometry::Point gravity(dim);
  for (int k = 0; k < dim; k++) gravity[k] = gravity_[k];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> flag(nfaces_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    AmanziGeometry::Point Kg = (*K_)[c] * gravity;
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);

      if (f < nfaces_owned && !flag[f]) {
        darcy_mass_flux[0][f] += rho * (Kg * normal);
        flag[f] = 1;
      }
    }
  }
}


/* ******************************************************************
* Updates global Darcy vector calculated by a discretization method.                                             
****************************************************************** */
void Matrix_MFD::AddGravityFluxes_DarcyFlux(Epetra_MultiVector& darcy_mass_flux,
                                         RelativePermeability& rel_perm) 
{
  double rho = *(S_->GetScalarData("fluid_density"));

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> flag(nfaces_wghost, 0);

  Epetra_MultiVector& Krel_cells = *rel_perm.Krel().ViewComponent("cell");
  Epetra_MultiVector& Krel_faces = *rel_perm.Krel().ViewComponent("face", true);
  int method = rel_perm.method();

  const Epetra_Vector& gravity_ = *S_->GetConstantVectorData("gravity");
  int dim = mesh_->space_dimension();

  AmanziGeometry::Point gravity(dim);
  for (int k = 0; k < dim; k++) gravity[k] = gravity_[k];

  for (int c = 0; c < ncells_owned; c++) {
    AmanziGeometry::Point Kg = (*K_)[c] * gravity;
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      std::vector<double>& krel = rel_perm.Krel_amanzi()[c];

      if (f < nfaces_owned && !flag[f]) {
        if (method == FLOW_RELATIVE_PERM_NONE) {
          darcy_mass_flux[0][f] += rho * (Kg * normal);
        } else if (method == FLOW_RELATIVE_PERM_CENTERED) {
          darcy_mass_flux[0][f] += rho * (Kg * normal) * Krel_cells[0][c];
        } else if (method == FLOW_RELATIVE_PERM_AMANZI) {
          darcy_mass_flux[0][f] += rho * (Kg * normal) * krel[n];
        } else {
          darcy_mass_flux[0][f] += rho * (Kg * normal) * Krel_faces[0][f];
        }
        flag[f] = 1;
      }
    }
  }
}




}  // namespace AmanziFlow
}  // namespace Amanzi

