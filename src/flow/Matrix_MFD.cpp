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
#include "Matrix_MFD.hpp"

namespace Amanzi {
namespace AmanziFlow {

Matrix_MFD::~Matrix_MFD()
{
  // if (MLprec->IsPreconditionerComputed()) {
  //  MLprec->DestroyPreconditioner();
  //  delete MLprec;
  // }
}


/* ******************************************************************
* Calculate elemental inverse mass matrices. 
* WARNING: The original Aff matrices are destroyed.                                            
****************************************************************** */
void Matrix_MFD::createMFDmassMatrices(int mfd3d_method, std::vector<WhetStone::Tensor>& K)
{
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Aff_cells.clear();
  for (int c = 0; c < K.size(); c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double> Bff(nfaces, nfaces);

    if (mfd3d_method == AmanziFlow::FLOW_MFD3D_HEXAHEDRA_MONOTONE) {
       if ((nfaces == 6 && dim == 3) || (nfaces == 4 && dim == 2))
         mfd.darcy_mass_inverse_hex(c, K[c], Bff);
       else
         mfd.darcy_mass_inverse(c, K[c], Bff);
         // mfd.darcy_mass_inverse_diagonal(c, K[c], Bff);
    } else {
       mfd.darcy_mass_inverse(c, K[c], Bff);
    }

    Aff_cells.push_back(Bff);
  }
}


/* ******************************************************************
* Calculate elemental stiffness matrices.                                            
****************************************************************** */
void Matrix_MFD::createMFDstiffnessMatrices(int mfd3d_method,
                                            std::vector<WhetStone::Tensor>& K,
                                            Epetra_Vector& Krel_faces)
{
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Aff_cells.clear();
  Afc_cells.clear();
  Acf_cells.clear();
  Acc_cells.clear();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double> Bff(nfaces, nfaces);
    Epetra_SerialDenseVector Bcf(nfaces), Bfc(nfaces);

    if (mfd3d_method == AmanziFlow::FLOW_MFD3D_HEXAHEDRA_MONOTONE) {
       if ((nfaces == 6 && dim == 3) || (nfaces == 4 && dim == 2))
         mfd.darcy_mass_inverse_hex(c, K[c], Bff);
       else
         mfd.darcy_mass_inverse(c, K[c], Bff);
         // mfd.darcy_mass_inverse_diagonal(c, K[c], Bff);
    } else {
       mfd.darcy_mass_inverse(c, K[c], Bff);
    }

    for (int n = 0; n < nfaces; n++)
      for (int m = 0; m < nfaces; m++) Bff(m, n) *= Krel_faces[faces[m]];

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

    Aff_cells.push_back(Bff);  // This the only place where memory can be allocated.
    Afc_cells.push_back(Bfc);
    Acf_cells.push_back(Bcf);
    Acc_cells.push_back(matsum);
  }
}


/* ******************************************************************
* May be used in the future.                                            
****************************************************************** */
void Matrix_MFD::rescaleMFDstiffnessMatrices(const Epetra_Vector& old_scale,
                                             const Epetra_Vector& new_scale)
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells; c++) {
    Teuchos::SerialDenseMatrix<int, double>& Bff = Aff_cells[c];
    Epetra_SerialDenseVector& Bcf = Acf_cells[c];

    int n = Bff.numRows();
    double scale = old_scale[c] / new_scale[c];

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) Bff(i, j) *= scale;
      Bcf(i) *= scale;
    }
    Acc_cells[c] *= scale;
  }
}


/* ******************************************************************
* Simply allocates memory.                                           
****************************************************************** */
void Matrix_MFD::createMFDrhsVectors()
{
  Ff_cells.clear();
  Fc_cells.clear();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Epetra_SerialDenseVector Ff(nfaces);  // Entries are initilaized to 0.0.
    double Fc = 0.0;

    Ff_cells.push_back(Ff);
    Fc_cells.push_back(Fc);
  }
}


/* ******************************************************************
* Applies boundary conditions to elemental stiffness matrices and
* creates elemental rigth-hand-sides.                                           
****************************************************************** */
void Matrix_MFD::applyBoundaryConditions(
    std::vector<int>& bc_markers, std::vector<double>& bc_values)
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double>& Bff = Aff_cells[c];  // B means elemental.
    Epetra_SerialDenseVector& Bfc = Afc_cells[c];
    Epetra_SerialDenseVector& Bcf = Acf_cells[c];

    Epetra_SerialDenseVector& Ff = Ff_cells[c];
    double& Fc = Fc_cells[c];

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      if (bc_markers[f] == FLOW_BC_FACE_PRESSURE ||
          bc_markers[f] == FLOW_BC_FACE_HEAD) {
        for (int m = 0; m < nfaces; m++) {
          Ff[m] -= Bff(m, n) * bc_values[f];
          Bff(n, m) = Bff(m, n) = 0.0;
        }
        Fc -= Bcf(n) * bc_values[f];
        Bcf(n) = Bfc(n) = 0.0;

        Bff(n, n) = 1.0;
        Ff[n] = bc_values[f];
      } else if (bc_markers[f] == FLOW_BC_FACE_FLUX) {
        Ff[n] -= bc_values[f] * mesh_->face_area(f);
      }
    }
  }
}


/* ******************************************************************
* Initialize Trilinos matrices. It must be called only once. 
* If matrix is non-symmetric, we generate transpose of the matrix 
* block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.   
****************************************************************** */
void Matrix_MFD::symbolicAssembleGlobalMatrices(const Epetra_Map& super_map)
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
  Acc = Teuchos::rcp(new Epetra_Vector(cmap));
  Acf = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cf_graph));
  Aff = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  Sff = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  Aff->GlobalAssemble();
  Sff->GlobalAssemble();

  if (flag_symmetry_)
    Afc = Acf;
  else
    Afc = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cf_graph));

  rhs = Teuchos::rcp(new Epetra_Vector(super_map));
  rhs_cells = Teuchos::rcp(FS->createCellView(*rhs));
  rhs_faces = Teuchos::rcp(FS->createFaceView(*rhs));
}


/* ******************************************************************
* Convert elemental mass matrices into stiffness matrices and 
* assemble them into four global matrices. 
* We need an auxiliary GHOST-based vector to assemble the RHS.
****************************************************************** */
void Matrix_MFD::assembleGlobalMatrices()
{
  Aff->PutScalar(0.0);

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
    (*Acc)[c] = Acc_cells[c];
    (*Acf).ReplaceMyValues(c, nfaces, Acf_cells[c].Values(), faces_LID);
    (*Aff).SumIntoGlobalValues(nfaces, faces_GID, Aff_cells[c].values());

    if (!flag_symmetry_)
        (*Afc).ReplaceMyValues(c, nfaces, Afc_cells[c].Values(), faces_LID);
  }
  (*Aff).GlobalAssemble();

  // We repeat some of the loops for code clarity.
  Epetra_Vector rhs_faces_wghost(fmap_wghost);

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    (*rhs_cells)[c] = Fc_cells[c];
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      rhs_faces_wghost[f] += Ff_cells[c][n];
    }
  }
  FS->combineGhostFace2MasterFace(rhs_faces_wghost, Add);

  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces; f++) (*rhs_faces)[f] = rhs_faces_wghost[f];
}


/* ******************************************************************
* Compute the face Schur complement of 2x2 block matrix.
****************************************************************** */
void Matrix_MFD::computeSchurComplement(
    std::vector<int>& bc_markers, std::vector<double>& bc_values)
{
  Sff->PutScalar(0.0);

  AmanziMesh::Entity_ID_List faces_LID;
  std::vector<int> dirs;
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces_LID, &dirs);
    int nfaces = faces_LID.size();
    Epetra_SerialDenseMatrix Schur(nfaces, nfaces);

    Epetra_SerialDenseVector& Bcf = Acf_cells[c];
    Epetra_SerialDenseVector& Bfc = Afc_cells[c];

    for (int n = 0; n < nfaces; n++) {
      for (int m = 0; m < nfaces; m++) {
        Schur(n, m) = Aff_cells[c](n, m) - Bfc[n] * Bcf[m] / (*Acc)[c];
      }
    }

    for (int n = 0; n < nfaces; n++) {  // Symbolic boundary conditions
      int f = faces_LID[n];
      if (bc_markers[f] == FLOW_BC_FACE_PRESSURE ||
          bc_markers[f] == FLOW_BC_FACE_HEAD) {
        for (int m = 0; m < nfaces; m++) Schur(n, m) = Schur(m, n) = 0.0;
        Schur(n, n) = 1.0;
      }
    }

    Epetra_IntSerialDenseVector faces_GID(nfaces);
    for (int n = 0; n < nfaces; n++) faces_GID[n] = (*Acf).ColMap().GID(faces_LID[n]);
    (*Sff).SumIntoGlobalValues(faces_GID, Schur);
  }
  (*Sff).GlobalAssemble();
}


/* ******************************************************************
* Linear algebra operations with matrices: r = f - A * x                                                 
****************************************************************** */
double Matrix_MFD::computeResidual(const Epetra_Vector& solution, Epetra_Vector& residual)
{
  Apply(solution, residual);
  residual.Update(1.0, *rhs, -1.0);

  double norm_residual;
  residual.Norm2(&norm_residual);
  return norm_residual;
}


/* ******************************************************************
* Linear algebra operations with matrices: r = A * x - f                                                 
****************************************************************** */
double Matrix_MFD::computeNegativeResidual(const Epetra_Vector& solution, Epetra_Vector& residual)
{
  Apply(solution, residual);
  residual.Update(-1.0, *rhs, 1.0);

  double norm_residual;
  residual.Norm2(&norm_residual);
  return norm_residual;
}


/* ******************************************************************
* Initialization of the preconditioner                                                 
****************************************************************** */
void Matrix_MFD::init_ML_preconditioner(Teuchos::ParameterList& ML_list_)
{
  ML_list = ML_list_;
  MLprec = new ML_Epetra::MultiLevelPreconditioner(*Sff, ML_list, false);
}


/* ******************************************************************
* Rebuild ML preconditioner.                                                 
****************************************************************** */
void Matrix_MFD::update_ML_preconditioner()
{
  if (MLprec->IsPreconditionerComputed()) MLprec->DestroyPreconditioner();
  MLprec->SetParameterList(ML_list);
  MLprec->ComputePreconditioner();
}


/* ******************************************************************
* Parallel matvec product Aff_cells * X.                                              
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
  ierr  = (*Aff).Multiply(false, Xf, Yf);
  ierr |= (*Afc).Multiply(true, Xc, Tf);  // Afc is kept in transpose form
  Yf.Update(1.0, Tf, 1.0);

  // Cell unknowns:  Yc = Acf * Xf + Acc * Xc
  ierr |= (*Acf).Multiply(false, Xf, Yc);  // It performs the required parallel communications.
  ierr |= Yc.Multiply(1.0, *Acc, Xc, 1.0);

  if (ierr) {
    Errors::Message msg("Matrix_MFD::Apply has failed to calculate Y = inv(A) * X.");
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
  ierr  = Tc.ReciprocalMultiply(1.0, *Acc, Xc, 0.0);
  ierr |= (*Afc).Multiply(true, Tc, Tf);  // Afc is kept in transpose form
  Tf.Update(1.0, Xf, -1.0);

  // Solve the Schur complement system Sff * Yf = Tf.
  MLprec->ApplyInverse(Tf, Yf);

  // BACKWARD SUBSTITUTION:  Yc = inv(Acc) (Xc - Acf Yf)
  ierr |= (*Acf).Multiply(false, Yf, Tc);  // It performs the required parallel communications.
  Tc.Update(1.0, Xc, -1.0);
  ierr |= Yc.ReciprocalMultiply(1.0, *Acc, Tc, 0.0);

  if (ierr) {
    Errors::Message msg("Matrix_MFD::ApplyInverse has failed in calculating y = A*x.");
    Exceptions::amanzi_throw(msg);
  }
  delete [] fvec_ptrs;
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
void Matrix_MFD::deriveDarcyMassFlux(const Epetra_Vector& solution,
                                     const Epetra_Import& face_importer,
                                     Epetra_Vector& darcy_mass_flux)
{
  Epetra_Vector* solution_faces = FS->createFaceView(solution);
#ifdef HAVE_MPI
  Epetra_Vector solution_faces_wghost(mesh_->face_map(true));
  solution_faces_wghost.Import(*solution_faces, face_importer, Insert);
#else
  Epetra_Vector& solution_faces_wghost = *solution_faces;
#endif

  AmanziMesh::Entity_ID_List faces;
  std::vector<double> dp;
  std::vector<int> dirs;

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  std::vector<int> flag(nfaces_wghost, 0);

  for (int c = 0; c < ncells; c++) {
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
        for (int m = 0; m < nfaces; m++) s += Aff_cells[c](n, m) * dp[m];
        darcy_mass_flux[f] = s * dirs[n];
        flag[f] = 1;
      }
    }
  }
}


/* ******************************************************************
* Derive Darcy velocity in cells. 
* WARNING: It cannot be consistent with the Darcy flux.                                                 
****************************************************************** */
void Matrix_MFD::deriveDarcyVelocity(const Epetra_Vector& darcy_flux,
                                     const Epetra_Import& face_importer,
                                     Epetra_MultiVector& darcy_velocity) const
{
#ifdef HAVE_MPI
  Epetra_Vector darcy_flux_wghost(mesh_->face_map(true));
  darcy_flux_wghost.Import(darcy_flux, face_importer, Insert);
#else
  Epetra_Vector& darcy_flux_wghost = darcy_flux;
#endif

  Teuchos::LAPACK<int, double> lapack;

  int dim = mesh_->space_dimension();
  Teuchos::SerialDenseMatrix<int, double> matrix(dim, dim);
  double rhs_cell[dim];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int i = 0; i < dim; i++) rhs_cell[i] = 0.0;
    matrix.putScalar(0.0);

    for (int n = 0; n < nfaces; n++) {  // populate least-square matrix
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double area = mesh_->face_area(f);

      for (int i = 0; i < dim; i++) {
        rhs_cell[i] += normal[i] * darcy_flux_wghost[f];
        matrix(i, i) += normal[i] * normal[i];
        for (int j = i+1; j < dim; j++) {
          matrix(j, i) = matrix(i, j) += normal[i] * normal[j];
        }
      }
    }

    int info;
    lapack.POSV('U', dim, 1, matrix.values(), dim, rhs_cell, dim, &info);

    for (int i = 0; i < dim; i++) darcy_velocity[i][c] = rhs_cell[i];
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

