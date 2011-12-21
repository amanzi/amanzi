/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include "Epetra_FECrsGraph.h"

#include "Flow_PK.hpp"
#include "Matrix_MFD.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Calculate elemental inverse mass matrices.                                            
****************************************************************** */
void Matrix_MFD::createMFDmassMatrices(std::vector<WhetStone::Tensor>& K)
{
  WhetStone::MFD3D mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;

  Minv_cells.clear();
  for (int c=0; c<K.size(); c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();
    Teuchos::SerialDenseMatrix<int, double> Wc(nfaces, nfaces);

    mfd.darcy_mass_inverse(c, K[c], Wc);
    Minv_cells.push_back(Wc);
  }
}


/* ******************************************************************
* Calculate elemental inverse mass matrices.                                            
****************************************************************** */
void Matrix_MFD::createMFDstiffnessMatrices(std::vector<WhetStone::Tensor>& K)
{
  WhetStone::MFD3D mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;

  Aff_cells.clear();
  Acc_cells.clear();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c<ncells; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double> Bff(nfaces, nfaces);
    Epetra_SerialDenseVector Bcf(nfaces);

    mfd.darcy_mass_inverse(c, K[c], Bff);

    double matsum = 0.0;  // elimination of mass matrix
    for (int n=0; n<nfaces; n++) {
      double colsum = 0.0;
      for (int m=0; m<nfaces; m++) colsum += Bff(n, m);
      Bcf(n) = -colsum;
      matsum += colsum;
    }

    Aff_cells.push_back(Bff);
    Acf_cells.push_back(Bcf);
    Acc_cells.push_back(matsum);
  }
}

/* ******************************************************************
* Applies boundary conditions to elemental stiffness matrices and
* creates elemental rigth-hand-sides.                                           
****************************************************************** */
void Matrix_MFD::applyBoundaryConditions(
    std::vector<int>& bc_markers, std::vector<double>& bc_values) 
{
  Ff_cells.clear();
  Fc_cells.clear();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;

  for (int c=0; c<ncells; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double>& Bff = Aff_cells[c];  // B means elemental.
    Epetra_SerialDenseVector& Bcf = Acf_cells[c];

    Epetra_SerialDenseVector Ff(nfaces);  // Entries are initilaized to 0.0.
    double Fc = 0.0;

    for (int n=0; n<nfaces; n++) {
      int f=faces[n];
      if (bc_markers[f] == FLOW_BC_FACE_PRESSURE || 
          bc_markers[f] == FLOW_BC_FACE_HEAD) {
        for (int m=0; m<nfaces; m++) { 
          Ff[n] -= Bff(n, m) * bc_values[f]; 
          Bff(n, m) = Bff(m, n) = 0.0;
        }
        Fc -= Bcf(n) * bc_values[f];
        Bcf(n) = 0.0;

        Bff(n, n) = 1.0;
        Ff(n) = bc_values[f]; 
      } else if (bc_markers[f] == FLOW_BC_FACE_FLUX) {
        // We assume for moment no-flux b.c.
      }
    }
    Ff_cells.push_back(Ff);
    Fc_cells.push_back(Fc);
  }
}


/* ******************************************************************
* Initialize Trilinos matrices. It must be called only once.                                           
****************************************************************** */
void Matrix_MFD::symbolicAssembleGlobalMatrices(const Epetra_Map& super_map)
{
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? FLOW_QUAD_FACES : FLOW_HEX_FACES;
  Epetra_CrsGraph cf_graph(Copy, cmap, super_map, avg_entries_row, true);
  Epetra_FECrsGraph ff_graph(Copy, fmap, 2*avg_entries_row);

  AmanziMesh::Entity_ID_List faces;
  int faces_LID[FLOW_MAX_FACES];  // Contigious memory is required.
  int faces_GID[FLOW_MAX_FACES];

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c<ncells; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n=0; n<nfaces; n++) {
      faces_LID[n] = faces[n];
      faces_GID[n] = fmap.GID(faces_LID[n]);
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

  //Aff->GlobalAssemble();
  //Sff->GlobalAssemble();
}


/* ******************************************************************
* Convert elemental mass matrices into stiffness matrices and 
* assemble them into four global matrices. 
* We need an auxiliary GHOST-based vector to assemble the RHS.
****************************************************************** */
void Matrix_MFD::assembleGlobalMatrices(Epetra_Vector& rhs)
{
  Aff->PutScalar(0.0);

  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(true);
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  AmanziMesh::Entity_ID_List faces;
  int faces_LID[FLOW_MAX_FACES];
  int faces_GID[FLOW_MAX_FACES];

  for (int c=0; c<ncells; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n=0; n<nfaces; n++) { 
      faces_LID[n] = faces[n];
      faces_GID[n] = fmap.GID(faces_LID[n]);
    }
    (*Acc)[c] = Acc_cells[c];
    (*Acf).ReplaceMyValues(c, nfaces, Acf_cells[c].Values(), faces_LID);
    (*Aff).SumIntoGlobalValues(nfaces, faces_GID, Aff_cells[c].values());
  }
  (*Aff).GlobalAssemble();

  // repeat some of the loops for clarity
#ifdef HAVE_MPI
  Epetra_Vector rhs_wghost(mesh_->face_map(true));
#else
  Epetra_Vector& rhs_wghost = rhs;
#endif
  for (int c=0; c<ncells; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    rhs[c] = Fc_cells[c];
    for (int n=0;  n<nfaces; n++) {
      int f = faces[n];
      rhs_wghost[f] += Ff_cells[c][n];
    }
  }
  FS->distribute_face_vector(rhs_wghost, Add);

#ifdef HAVE_MPI
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f<nfaces_owned; f++) rhs[ncells + f] = rhs_wghost[f];
#endif
}


/* ******************************************************************
* Compute the face Schur complement of 2x2 block matrix.
* Special structure of matrix Acf is exploited in this algorithm.                               
****************************************************************** */
void Matrix_MFD::computeSchurComplement()
{
  // Simple *Sff = *Aff does other stuff we do not want?
  int *offsets, *indices;
  double *source, *target;
  (*Aff).ExtractCrsDataPointers(offsets, indices, source);
  (*Sff).ExtractCrsDataPointers(offsets, indices, target);
  for (int n=0; n<(*Aff).NumMyNonzeros(); n++) target[n] = source[n];

  double *Bcf;
  int nfaces, *faces_LID;
  int ncells = (*Acc).Map().NumMyElements();

  for (int c=0; c<ncells; c++) {
    (*Acf).ExtractMyRowView(c, nfaces, Bcf, faces_LID);
    Epetra_SerialDenseMatrix Schur(nfaces, nfaces);

    for (int n=0; n<nfaces; n++) {
      for (int m=0; m<nfaces; m++) {
        Schur(n, m) = -Bcf[n] * Bcf[m] / (*Acc)[c];
      }
    }
    Epetra_IntSerialDenseVector faces_GID(nfaces);
    for (int n=0; n<nfaces; n++) faces_GID[n] = (*Acf).ColMap().GID(faces_LID[n]);
    (*Sff).SumIntoGlobalValues(faces_GID, Schur);
  }
  (*Sff).GlobalAssemble();  // Do we really need it? (lipnikov@lanl.gov)
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
  double **fvec_ptrs = new double*[X.NumVectors()];
  for (int i=0; i<nvectors; i++) fvec_ptrs[i] = cvec_ptrs[i] + ncells;

  Epetra_MultiVector Xc(View, cmap, cvec_ptrs, nvectors);
  Epetra_MultiVector Xf(View, fmap, fvec_ptrs, nvectors);

  // Create views Yc and Yf into the cell and face segments of Y.
  cvec_ptrs = Y.Pointers();
  for (int i=0; i<nvectors; i++) fvec_ptrs[i] = cvec_ptrs[i] + ncells;

  Epetra_MultiVector Yc(View, cmap, cvec_ptrs, nvectors);
  Epetra_MultiVector Yf(View, fmap, fvec_ptrs, nvectors);

  // Temporary cell and face vectors.
  Epetra_MultiVector Tc(cmap, nvectors);
  Epetra_MultiVector Tf(fmap, nvectors);

  // FORWARD ELIMINATION:  Tf = Xf - trans(Acf) inv(Acc) Xc
  Tc.ReciprocalMultiply(1.0, *Acc, Xc, 0.0);
  (*Acf).Multiply(true, Tc, Tf);  // this should do the required parallel comm
  Tf.Update(1.0, Xf, -1.0);

  // Solve the Schur complement system for Yf with Tf as the rhs using ML
  MLprec->ApplyInverse(Tf, Yf);

  // BACKWARD SUBSTITUTION: Yc = inv(Acc) (Xc - Acf Yf)
  Tf = Yf;
  (*Acf).Multiply(false, Tf, Tc);  // this should do the required parallel comm
  Tc.Update(1.0, Xc, -1.0);
  Yc.ReciprocalMultiply(1.0, *Acc, Tc, 0.0);

  delete [] cvec_ptrs, fvec_ptrs;
  return 0;
}


/* ******************************************************************
* Derive Darcy flux.                                                 
****************************************************************** */
void Matrix_MFD::deriveDarcyFlux(const Epetra_Vector& solution, 
                                 const Epetra_Vector& rhs, 
                                 const Epetra_Import& face_importer, 
                                 Epetra_Vector& darcy_flux)
{
  Epetra_Vector* solution_faces = FS->create_face_view(solution);
#ifdef HAVE_MPI
  Epetra_Vector solution_faces_wghost(mesh_->face_map(true));
  solution_faces_wghost.Import(*solution_faces, face_importer, Insert);
#else
  Epetra_Vector& solution_faces_wghost = solution_faces;
#endif

  AmanziMesh::Entity_ID_List faces;
  std::vector<double> dp;
  std::vector<int> dirs;

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c<ncells; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    dp.resize(nfaces);
    for (int n=0; n<nfaces; n++) { 
      int f = faces[n];
      dp[n] = rhs[c] + solution[c] - solution_faces_wghost[f];
    }

    mesh_->cell_get_face_dirs(c, &dirs);
    for (int n=0; n<nfaces; n++) {
      double s = 0.0;
      for (int m=0; m<nfaces; m++) s += Aff_cells[c](n, m) * dp[m];
      int f = faces[n];
      darcy_flux[f] = s * dirs[n];
    }
  }
}


/* ******************************************************************
* Derive Darcy velocity in cells. 
* WARNING: It cannot be consistent with the Darcy flux.                                                 
****************************************************************** */
void Matrix_MFD::deriveDarcyVelocity(const Epetra_Vector& darcy_flux, Epetra_MultiVector& darcy_velocity) const
{
  Teuchos::LAPACK<int, double> lapack;
  double rhs[FLOW_MAX_FACES];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int dim = mesh_->space_dimension();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c<ncells; c++) {
    mesh_->cell_get_faces(c, &faces); 
    mesh_->cell_get_face_dirs(c, &dirs);
    int nfaces = faces.size(); 

    Teuchos::SerialDenseMatrix<int, double> matrix(nfaces, nfaces); 

    for (int n=0; n<ncells; n++) {  // populate least-square matrix
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double area = mesh_->face_area(f);

      for (int i=0; i<dim; i++) {
        rhs[i] += normal[i] * darcy_flux[f] * dirs[f];
        matrix(i,i) += normal[i] * normal[i]; 
        for (int j=i+1; j<dim; j++) {
          matrix(j,i) = matrix(i,j) += normal[i] * normal[j];
        }
      }
    }

    int info;
    lapack.POSV('U', dim, 1, matrix.values(), dim, rhs, dim, &info); 

    for (int i=0; i<dim; i++) darcy_velocity[i][c] = rhs[i];
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

