/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <vector>

#include "Epetra_FECrsGraph.h"

#include "Flow_PK.hpp"
#include "Flow_constants.hpp"
#include "Matrix_MFD.hpp"


namespace Amanzi {
namespace AmanziFlow {

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
void Matrix_MFD::CreateMFDmassMatrices(int mfd3d_method, std::vector<WhetStone::Tensor>& K)
{
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Mff_cells_.clear();

  int ok;
  nokay_ = npassed_ = 0;

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double> Mff(nfaces, nfaces);

    if (mfd3d_method == AmanziFlow::FLOW_MFD3D_HEXAHEDRA_MONOTONE) {
      if ((nfaces == 6 && dim == 3) || (nfaces == 4 && dim == 2))
        ok = mfd.DarcyMassInverseHex(c, K[c], Mff);
      else
        ok = mfd.DarcyMassInverse(c, K[c], Mff);
    } else if (mfd3d_method == AmanziFlow::FLOW_MFD3D_TWO_POINT_FLUX) {
      ok = mfd.DarcyMassInverseDiagonal(c, K[c], Mff);
    } else if (mfd3d_method == AmanziFlow::FLOW_MFD3D_SUPPORT_OPERATOR) {
      ok = mfd.DarcyMassInverseSO(c, K[c], Mff);
    } else if (mfd3d_method == AmanziFlow::FLOW_MFD3D_OPTIMIZED) {
      ok = mfd.DarcyMassInverseOptimized(c, K[c], Mff);
    } else {
      ok = mfd.DarcyMassInverse(c, K[c], Mff);
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
* Calculate elemental stiffness matrices.                                            
****************************************************************** */
void Matrix_MFD::CreateMFDstiffnessMatrices(Epetra_Vector& Krel_cells,
                                            Epetra_Vector& Krel_faces,
                                            int method)
{
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Aff_cells_.clear();
  Afc_cells_.clear();
  Acf_cells_.clear();
  Acc_cells_.clear();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double>& Mff = Mff_cells_[c];
    Teuchos::SerialDenseMatrix<int, double> Bff(nfaces, nfaces);
    Epetra_SerialDenseVector Bcf(nfaces), Bfc(nfaces);

    if (method == FLOW_RELATIVE_PERM_NONE) {
      double* braw = Bff.values();
      double* mraw = Mff.values();
      for (int n = 0; n < nfaces * nfaces; n++) braw[n] = mraw[n];
    } else if (method == FLOW_RELATIVE_PERM_CENTERED ||
               method == FLOW_RELATIVE_PERM_EXPERIMENTAL) {  // centered permeability for diffusion
      for (int n = 0; n < nfaces; n++)
        for (int m = 0; m < nfaces; m++) Bff(m, n) = Mff(m, n) * Krel_cells[c];
    } else {
      for (int n = 0; n < nfaces; n++)
        for (int m = 0; m < nfaces; m++) Bff(m, n) = Mff(m, n) * Krel_faces[faces[m]];
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
void Matrix_MFD::CreateMFDrhsVectors()
{
  Ff_cells_.clear();
  Fc_cells_.clear();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Epetra_SerialDenseVector Ff(nfaces);  // Entries are initilaized to 0.0.
    double Fc = 0.0;

    Ff_cells_.push_back(Ff);
    Fc_cells_.push_back(Fc);
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
  std::vector<int> dirs;

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double>& Bff = Aff_cells_[c];  // B means elemental.
    Epetra_SerialDenseVector& Bfc = Afc_cells_[c];
    Epetra_SerialDenseVector& Bcf = Acf_cells_[c];

    Epetra_SerialDenseVector& Ff = Ff_cells_[c];
    double& Fc = Fc_cells_[c];

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      double value = bc_values[f][0];

      if (bc_model[f] == FLOW_BC_FACE_PRESSURE ||
          bc_model[f] == FLOW_BC_FACE_PRESSURE_SEEPAGE) {
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
      } else if (bc_model[f] == FLOW_BC_FACE_MIXED) {  // Not used yet
        double area = mesh_->face_area(f);
        Ff[n] += value * area;
        Bff(n, n) += bc_values[f][1] * area;
      }

      // Additional work required for seepage boundary condition.
      if (bc_model[f] == FLOW_BC_FACE_PRESSURE_SEEPAGE) {
        Fc -= bc_values[f][1] * mesh_->face_area(f);
      }
    }
  }
}


/* ******************************************************************
* Initialize Trilinos matrices. It must be called only once. 
* If matrix is non-symmetric, we generate transpose of the matrix 
* block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.   
****************************************************************** */
void Matrix_MFD::SymbolicAssembleGlobalMatrices(const Epetra_Map& super_map)
{
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? FLOW_QUAD_FACES : FLOW_HEX_FACES;
  Epetra_CrsGraph cf_graph(Copy, cmap, fmap_wghost, avg_entries_row, false);  // FIX (lipnikov@lanl.gov)
  Epetra_FECrsGraph ff_graph(Copy, fmap, 2*avg_entries_row);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int faces_LID[FLOW_MAX_FACES];  // Contigious memory is required.
  int faces_GID[FLOW_MAX_FACES];

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      faces_LID[n] = faces[n];
      faces_GID[n] = fmap_wghost.GID(faces_LID[n]);
    }
    cf_graph.InsertMyIndices(c, nfaces, faces_LID);
    ff_graph.InsertGlobalIndices(nfaces, faces_GID, nfaces, faces_GID);
  }
  cf_graph.FillComplete(fmap, cmap);
  ff_graph.GlobalAssemble();  // Symbolic graph is complete.

  // create global matrices
  Acc_ = Teuchos::rcp(new Epetra_Vector(cmap));
  Acf_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cf_graph));
  Aff_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  Sff_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  Aff_->GlobalAssemble();
  Sff_->GlobalAssemble();

  if (flag_symmetry_)
    Afc_ = Acf_;
  else
    Afc_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cf_graph));

  rhs_ = Teuchos::rcp(new Epetra_Vector(super_map));
  rhs_cells_ = Teuchos::rcp(FS->CreateCellView(*rhs_));
  rhs_faces_ = Teuchos::rcp(FS->CreateFaceView(*rhs_));
}


/* ******************************************************************
* Assebmle elemental mass matrices into four global matrices. 
* We need an auxiliary GHOST-based vector to assemble the RHS.
****************************************************************** */
void Matrix_MFD::AssembleGlobalMatrices()
{
  Aff_->PutScalar(0.0);

  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int faces_LID[FLOW_MAX_FACES];
  int faces_GID[FLOW_MAX_FACES];

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      faces_LID[n] = faces[n];
      faces_GID[n] = fmap_wghost.GID(faces_LID[n]);
    }
    (*Acc_)[c] = Acc_cells_[c];
    Acf_->ReplaceMyValues(c, nfaces, Acf_cells_[c].Values(), faces_LID);
    Aff_->SumIntoGlobalValues(nfaces, faces_GID, Aff_cells_[c].values());

    if (!flag_symmetry_)
        Afc_->ReplaceMyValues(c, nfaces, Afc_cells_[c].Values(), faces_LID);
  }
  Aff_->GlobalAssemble();

  // We repeat some of the loops for code clarity.
  Epetra_Vector rhs_faces_wghost(fmap_wghost);

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    (*rhs_cells_)[c] = Fc_cells_[c];
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      rhs_faces_wghost[f] += Ff_cells_[c][n];
    }
  }
  FS->CombineGhostFace2MasterFace(rhs_faces_wghost, Add);

  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces; f++) (*rhs_faces_)[f] = rhs_faces_wghost[f];
}


/* ******************************************************************
* Compute the face Schur complement of 2x2 block matrix.
****************************************************************** */
void Matrix_MFD::ComputeSchurComplement(
    std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values)
{
  Sff_->PutScalar(0.0);

  AmanziMesh::Entity_ID_List faces_LID;
  std::vector<int> dirs;
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces_LID, &dirs);
    int nfaces = faces_LID.size();
    Epetra_SerialDenseMatrix Schur(nfaces, nfaces);

    Epetra_SerialDenseVector& Bcf = Acf_cells_[c];
    Epetra_SerialDenseVector& Bfc = Afc_cells_[c];

    for (int n = 0; n < nfaces; n++) {
      for (int m = 0; m < nfaces; m++) {
        Schur(n, m) = Aff_cells_[c](n, m) - Bfc[n] * Bcf[m] / (*Acc_)[c];
      }
    }

    for (int n = 0; n < nfaces; n++) {  // Symbolic boundary conditions
      int f = faces_LID[n];
      if (bc_model[f] == FLOW_BC_FACE_PRESSURE ||
          bc_model[f] == FLOW_BC_FACE_PRESSURE_SEEPAGE) {
        for (int m = 0; m < nfaces; m++) Schur(n, m) = Schur(m, n) = 0.0;
        Schur(n, n) = 1.0;
      }
    }

    Epetra_IntSerialDenseVector faces_GID(nfaces);
    for (int n = 0; n < nfaces; n++) faces_GID[n] = (*Acf_).ColMap().GID(faces_LID[n]);
    (*Sff_).SumIntoGlobalValues(faces_GID, Schur);
  }
  (*Sff_).GlobalAssemble();
}


/* ******************************************************************
* Linear algebra operations with matrices: r = f - A * x                                                 
****************************************************************** */
double Matrix_MFD::ComputeResidual(const Epetra_Vector& solution, Epetra_Vector& residual)
{
  Apply(solution, residual);
  residual.Update(1.0, *rhs_, -1.0);

  double norm_residual;
  residual.Norm2(&norm_residual);
  return norm_residual;
}


/* ******************************************************************
* Linear algebra operations with matrices: r = A * x - f                                                 
****************************************************************** */
double Matrix_MFD::ComputeNegativeResidual(const Epetra_Vector& solution, Epetra_Vector& residual)
{
  Apply(solution, residual);
  residual.Update(-1.0, *rhs_, 1.0);

  double norm_residual;
  residual.Norm2(&norm_residual);
  return norm_residual;
}


/* ******************************************************************
* Initialization of the preconditioner                                                 
****************************************************************** */
void Matrix_MFD::InitPreconditioner(int method, Teuchos::ParameterList& prec_list)
{
  method_ = method;

  if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
    ML_list = prec_list;
    MLprec = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*Sff_, ML_list, false));
  } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) {
#ifdef HAVE_HYPRE
    hypre_ncycles = prec_list.get<int>("cycle applications", 5);  // Boomer AMG parameters
    hypre_nsmooth = prec_list.get<int>("smoother sweeps", 3);
    hypre_tol = prec_list.get<double>("tolerance", 0.0);
    hypre_strong_threshold = prec_list.get<double>("strong threshold", 0.0);
#endif
  } else if (method_ == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    ifp_plist_ = prec_list;
  }
}


/* ******************************************************************
* Rebuild the preconditioner.                                                 
****************************************************************** */
void Matrix_MFD::UpdatePreconditioner()
{
  if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
    if (MLprec->IsPreconditionerComputed()) MLprec->DestroyPreconditioner();
    MLprec->SetParameterList(ML_list);
    MLprec->ComputePreconditioner();
  } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) {
#ifdef HAVE_HYPRE
    IfpHypre_Sff_ = Teuchos::rcp(new Ifpack_Hypre(&*Sff_));

    Teuchos::RCP<FunctionParameter> functs[8];
    functs[0] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCoarsenType, 0));
    functs[1] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetPrintLevel, 0)); 
    functs[2] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetNumSweeps, hypre_nsmooth));
    functs[3] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxIter, hypre_ncycles));
    functs[4] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetRelaxType, 6)); 
    functs[5] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetStrongThreshold, hypre_strong_threshold)); 
    functs[6] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetTol, hypre_tol)); 
    functs[7] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCycleType, 1));  

    Teuchos::ParameterList hypre_list("Preconditioner List");
    // hypre_list.set("Solver", PCG);
    hypre_list.set("Preconditioner", BoomerAMG);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 8);
    hypre_list.set<Teuchos::RCP<FunctionParameter>*>("Functions", functs); 

    IfpHypre_Sff_->SetParameters(hypre_list);
    IfpHypre_Sff_->Initialize();
    IfpHypre_Sff_->Compute();
#endif
  } else if (method_ == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    Ifpack factory;
    std::string prectype("ILU");
    int ovl = ifp_plist_.get<int>("overlap", 0);
    ifp_plist_.set<std::string>("schwarz: combine mode", "Add");
    ifp_prec_ = Teuchos::rcp(factory.Create(prectype, &*Sff_, ovl));
    ifp_prec_->SetParameters(ifp_plist_);
    ifp_prec_->Initialize();
    ifp_prec_->Compute();
  }
}


/* ******************************************************************
* Parallel matvec product A * X.                                              
****************************************************************** */
int Matrix_MFD::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nvectors = X.NumVectors();

  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);

  // Create views Xc and Xf into the cell and face segments of X.
  double **cvec_ptrs = X.Pointers();
  double **fvec_ptrs = new double*[nvectors];
  for (int i = 0; i < nvectors; i++) fvec_ptrs[i] = cvec_ptrs[i] + ncells;

  Epetra_MultiVector Xc(View, cmap, cvec_ptrs, nvectors);
  Epetra_MultiVector Xf(View, fmap, fvec_ptrs, nvectors);

  // Create views Yc and Yf into the cell and face segments of Y.
  cvec_ptrs = Y.Pointers();
  for (int i = 0; i < nvectors; i++) fvec_ptrs[i] = cvec_ptrs[i] + ncells;

  Epetra_MultiVector Yc(View, cmap, cvec_ptrs, nvectors);
  Epetra_MultiVector Yf(View, fmap, fvec_ptrs, nvectors);

  // Face unknowns:  Yf = Aff * Xf + Afc * Xc
  int ierr;
  Epetra_MultiVector Tf(fmap, nvectors);
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
  delete [] fvec_ptrs;
  return 0;
}


/* ******************************************************************
* The OWNED cell-based and face-based d.o.f. are packed together into 
* the X and Y Epetra vectors, with the cell-based in the first part.
*
* WARNING: When invoked by AztecOO the arguments X and Y may be 
* aliased: possibly the same object or different views of the same 
* underlying data. Thus, we do not assign to Y until the end.                                              
****************************************************************** */
int Matrix_MFD::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nvectors = X.NumVectors();

  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);

  // Create views Xc and Xf into the cell and face segments of X.
  double **cvec_ptrs = X.Pointers();
  double **fvec_ptrs = new double*[nvectors];
  for (int i = 0; i < nvectors; i++) fvec_ptrs[i] = cvec_ptrs[i] + ncells;

  Epetra_MultiVector Xc(View, cmap, cvec_ptrs, nvectors);
  Epetra_MultiVector Xf(View, fmap, fvec_ptrs, nvectors);

  // Create views Yc and Yf into the cell and face segments of Y.
  cvec_ptrs = Y.Pointers();
  for (int i = 0; i < nvectors; i++) fvec_ptrs[i] = cvec_ptrs[i] + ncells;

  Epetra_MultiVector Yc(View, cmap, cvec_ptrs, nvectors);
  Epetra_MultiVector Yf(View, fmap, fvec_ptrs, nvectors);

  // Temporary cell and face vectors.
  Epetra_MultiVector Tc(cmap, nvectors);
  Epetra_MultiVector Tf(fmap, nvectors);

  // FORWARD ELIMINATION:  Tf = Xf - Afc inv(Acc) Xc
  int ierr;
  ierr  = Tc.ReciprocalMultiply(1.0, *Acc_, Xc, 0.0);
  ierr |= (*Afc_).Multiply(true, Tc, Tf);  // Afc is kept in transpose form
  Tf.Update(1.0, Xf, -1.0);

  // Solve the Schur complement system Sff * Yf = Tf.
  if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
    MLprec->ApplyInverse(Tf, Yf);
  } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) { 
#ifdef HAVE_HYPRE
    ierr |= IfpHypre_Sff_->ApplyInverse(Tf, Yf);
#endif
  } else if (method_ == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    ifp_prec_->ApplyInverse(Tf, Yf);
  }

  // BACKWARD SUBSTITUTION:  Yc = inv(Acc) (Xc - Acf Yf)
  ierr |= (*Acf_).Multiply(false, Yf, Tc);  // It performs the required parallel communications.
  Tc.Update(1.0, Xc, -1.0);
  ierr |= Yc.ReciprocalMultiply(1.0, *Acc_, Tc, 0.0);

  if (ierr) {
    Errors::Message msg("Matrix_MFD::ApplyInverse has failed in calculating y = inv(A)*x.");
    Exceptions::amanzi_throw(msg);
  }
  delete [] fvec_ptrs;
  return 0;
}


/* ******************************************************************
* Reduce the pressure-lambda-system to lambda-system via ellimination
* of the known pressure. Structure of the global system is preserved
* but off-diagola blocks are zeroed-out.                                               
****************************************************************** */
int Matrix_MFD::ReduceGlobalSystem2LambdaSystem(Epetra_Vector& u)
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);

  // Create views of u and rhs
  double **c_ptrs = u.Pointers();
  double **f_ptrs = rhs_->Pointers();

  Epetra_Vector uc(View, cmap, c_ptrs[0]);
  Epetra_Vector gf(View, fmap, f_ptrs[0] + ncells);

  // Update RHS: rhs = rhs - Afc * uc
  Epetra_Vector tf(fmap);
  Afc_->Multiply(true, uc, tf);  // Afc is kept in the transpose form.
  gf.Update(-1.0, tf, 1.0);

  // Decomple pressure-lambda system
  *Sff_ = *Aff_;
  Afc_->PutScalar(0.0);
  Acf_->PutScalar(0.0);

  return 0;
}


/* ******************************************************************
* WARNING: Routines requires original mass matrices (Aff_cells), i.e.
* before boundary conditions were imposed.
*
* WARNING: Since diffusive flux is not continuous, we derive it only
* once (using flag) and in exactly the same manner as in routine
* Flow_PK::addGravityFluxes_DarcyFlux.
****************************************************************** */
void Matrix_MFD::DeriveDarcyMassFlux(const Epetra_Vector& solution,
                                     const Epetra_Import& face_importer,
                                     Epetra_Vector& darcy_mass_flux)
{
  Teuchos::RCP<Epetra_Vector> solution_faces = Teuchos::rcp(FS->CreateFaceView(solution));
#ifdef HAVE_MPI
  Epetra_Vector solution_faces_wghost(mesh_->face_map(true));
  solution_faces_wghost.Import(*solution_faces, face_importer, Insert);
#else
  Epetra_Vector& solution_faces_wghost = *solution_faces;
#endif

  AmanziMesh::Entity_ID_List faces;
  std::vector<double> dp;
  std::vector<int> dirs;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  std::vector<int> flag(nfaces_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    dp.resize(nfaces);
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      dp[n] = solution[c] - solution_faces_wghost[f];
    }

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      if (f < nfaces_owned && !flag[f]) {
        double s(0.0);
        for (int m = 0; m < nfaces; m++) s += Aff_cells_[c](n, m) * dp[m];
        darcy_mass_flux[f] = s * dirs[n];
        flag[f] = 1;
      }
    }
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

