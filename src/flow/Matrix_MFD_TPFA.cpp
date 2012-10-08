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

#include "Flow_constants.hpp"
#include "Matrix_MFD_TPFA.hpp"

namespace Amanzi {
namespace AmanziFlow {

template<class T>
int FindPosition(const std::vector<T>& v, const T& value) {
  for (int i = 0; i < v.size(); i++) 
     if (v[i] == value) return i;
  return -1;
}


/* ******************************************************************
* Calculate elemental stiffness matrices.                                            
****************************************************************** */
void Matrix_MFD_TPFA::CreateMFDstiffnessMatrices(Epetra_Vector& Krel_cells,
                                                 Epetra_Vector& Krel_faces)
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

    if (Krel_cells[c] == 1.0) {
      for (int n = 0; n < nfaces; n++) Bff(n, n) = Mff(n, n) * Krel_faces[faces[n]];
    } else {
      for (int n = 0; n < nfaces; n++) Bff(n, n) = Mff(n, n) * Krel_cells[c];
    }

    double matsum = 0.0;  // elimination of mass matrix
    for (int n = 0; n < nfaces; n++) {
      double rowsum = Bff(n, n), colsum = Bff(n, n);
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
* Initialize Trilinos matrices. It must be called only once. 
* If matrix is non-symmetric, we generate transpose of the matrix 
* block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.   
****************************************************************** */
void Matrix_MFD_TPFA::SymbolicAssembleGlobalMatrices(const Epetra_Map& super_map)
{ 
  Matrix_MFD::SymbolicAssembleGlobalMatrices(super_map);

  const Epetra_Map& fmap = mesh_->face_map(false);
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? FLOW_QUAD_FACES : FLOW_HEX_FACES;
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
  Dff_ = Teuchos::rcp(new Epetra_Vector(fmap));
  Spp_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, pp_graph));
  Spp_->GlobalAssemble();
}


/* ******************************************************************
* Convert elemental mass matrices into stiffness matrices and 
* assemble them into four global matrices. 
* We need an auxiliary GHOST-based vector to assemble the RHS.
****************************************************************** */
void Matrix_MFD_TPFA::AssembleGlobalMatrices()
{
  Matrix_MFD::AssembleGlobalMatrices();

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  Dff_->PutScalar(0.0);
  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      (*Dff_)[f] += Aff_cells_[c](n, n);
    }
  }

  // convert right-hand side to a cell-based vector
  const Epetra_Map& cmap = mesh_->cell_map(false);
  Epetra_Vector Tc(cmap);
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  for (int f = 0; f < nfaces; f++) (*rhs_faces_)[f] /= (*Dff_)[f];
  (*Acf_).Multiply(false, *rhs_faces_, Tc);
  for (int c = 0; c < ncells; c++) (*rhs_cells_)[c] -= Tc[c];
  
  rhs_faces_->PutScalar(0.0);

  // create a ghost version of Acc
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  Epetra_Vector Dcc(cmap_wghost);

  for (int c = 0; c < ncells; c++) Dcc[c] = (*Acc_)[c];
  FS->CopyMasterFace2GhostFace(Dcc);

  AmanziMesh::Entity_ID_List cells;
  int cells_GID[2];
  double Acf_copy[2];
  
  Spp_->PutScalar(0.0);
  for (int f = 0; f < nfaces; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int mcells = cells.size();

    // populate face-based matrix.
    Teuchos::SerialDenseMatrix<int, double> Bpp(mcells, mcells);
    for (int n = 0; n < mcells; n++) {
      int c = cells[n];
      cells_GID[n] = cmap_wghost.GID(c);

      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int i = FindPosition<int>(faces, f);
      Acf_copy[n] = Acf_cells_[c][i];
      Bpp(n, n) = Dcc[c] / faces.size();
    }

    for (int n = 0; n < mcells; n++) {
      for (int m = n; m < mcells; m++) {
        Bpp(n, m) -= Acf_copy[n] * Acf_copy[m] / (*Dff_)[f];
        Bpp(m, n) = Bpp(n, m);
      }        
    }
    
    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Bpp.values());
  }
  (*Spp_).GlobalAssemble();
}


/* ******************************************************************
* Compute the face Schur complement of 2x2 block matrix.
****************************************************************** */
void Matrix_MFD_TPFA::ComputeSchurComplement(
    std::vector<int>& bc_markers, std::vector<double>& bc_values)
{
}


/* ******************************************************************
* Initialization of the preconditioner                                                 
****************************************************************** */
void Matrix_MFD_TPFA::InitPreconditioner(int method, Teuchos::ParameterList& prec_list)
{
  method_ = method;

  if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
    ML_list = prec_list;
    MLprec = new ML_Epetra::MultiLevelPreconditioner(*Spp_, ML_list, false);
  } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) {
#ifdef HAVE_HYPRE
    // read some boomer amg parameters
    hypre_ncycles = prec_list.get<int>("cycle applications", 5);
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
void Matrix_MFD_TPFA::UpdatePreconditioner()
{
  if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
    if (MLprec->IsPreconditionerComputed()) MLprec->DestroyPreconditioner();
    MLprec->SetParameterList(ML_list);
    MLprec->ComputePreconditioner();
  } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) {
#ifdef HAVE_HYPRE
    IfpHypre_Spp_ = Teuchos::rcp(new Ifpack_Hypre(&*Spp_));
    Teuchos::RCP<FunctionParameter> functs[8];
    functs[0] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCoarsenType, 0));
    functs[1] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetPrintLevel, 0)); 
    functs[2] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetNumSweeps, hypre_nsmooth));
    functs[3] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxIter, hypre_ncycles));
    functs[4] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetRelaxType, 6)); 
    functs[5] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetStrongThreshold, hypre_strong_threshold)); 
    functs[6] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetTol, hypre_tol)); 
    functs[7] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCycleType, 1));  

    Teuchos::ParameterList hypre_list;
    // hypre_list.set("Solver", PCG);
    hypre_list.set("Preconditioner", BoomerAMG);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 8);
    hypre_list.set<Teuchos::RCP<FunctionParameter>*>("Functions", functs); 

    IfpHypre_Spp_->SetParameters(hypre_list);
    IfpHypre_Spp_->Initialize();
    IfpHypre_Spp_->Compute();
#endif
  } else if (method_ == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    Ifpack factory;
    std::string prectype("ILU");
    int ovl = ifp_plist_.get<int>("overlap", 0);
    ifp_plist_.set<std::string>("schwarz: combine mode", "Add");
    ifp_prec_ = Teuchos::rcp(factory.Create(prectype, &*Spp_, ovl));
    ifp_prec_->SetParameters(ifp_plist_);
    ifp_prec_->Initialize();
    ifp_prec_->Compute();
  }
}


/* ******************************************************************
* Parallel matvec product Spp * Xc.                                              
****************************************************************** */
int Matrix_MFD_TPFA::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nvectors = X.NumVectors();

  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);

  // Create views Xc into the cell segments of X.
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

  int ierr = (*Spp_).Multiply(false, Xc, Yc);
 
  if (ierr) {
    Errors::Message msg("Matrix_MFD_TPFA::Apply has failed to calculate y = A*x.");
    Exceptions::amanzi_throw(msg);
  }

  // Yf = Xf;
  Yf.PutScalar(0.0);

  delete [] fvec_ptrs;
  return 0;
}


/* ******************************************************************
* The OWNED cell-based and face-based d.o.f. are packed together into 
* the X and Y Epetra vectors, with the cell-based in the first part.                                           
****************************************************************** */
int Matrix_MFD_TPFA::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nvectors = X.NumVectors();

  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);

  // Create views Xc into the cell segments of X.
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

  // Solve the Schur complement system Spp * Yc = Xc. Since AztecOO may
  // use the same memory for X and Y, we introduce auxiliaty vector Tc.
  int ierr = 0;
  Epetra_Vector Tc(cmap);
  if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
    MLprec->ApplyInverse(Xc, Tc);
  } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) { 
#ifdef HAVE_HYPRE
    ierr = IfpHypre_Spp_->ApplyInverse(Xc, Tc);
#endif
  } else if (method_ == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    ifp_prec_->ApplyInverse(Xc, Tc);
  }
  Yc = Tc;

  if (ierr) {
    Errors::Message msg("Matrix_MFD_TPFA::ApplyInverse has failed in calculating y = inv(A)*x.");
    Exceptions::amanzi_throw(msg);
  }

  Yf = Xf;

  delete [] fvec_ptrs;
  return 0;
}

}  // namespace AmanziFlow
}  // namespace Amanzi


