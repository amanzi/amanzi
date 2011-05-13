#include "DiffusionMatrix.hpp"

#include "Epetra_FECrsGraph.h"
#include "MimeticHexLocal.hpp"

#include <iostream>

using namespace Amanzi;
using namespace AmanziMesh;

DiffusionMatrix::DiffusionMatrix(const Teuchos::RCP<Amanzi::AmanziMesh::Mesh> &mesh,
                                 const std::vector<int> &dir_faces) : mesh_(mesh), dir_faces_(dir_faces)
{
  const Epetra_Map& cell_map = CellMap(false); // owned cells
  const Epetra_Map& face_map = FaceMap(false); // owned faces
  const Epetra_Map& face_map_use = FaceMap(true);  // all used faces

  // Create graphs for the cell-face and face-face matrices.
  int l_indices[6], g_indices[6];
  Epetra_CrsGraph cf_graph(Copy, cell_map, face_map_use, 6, true);
  Epetra_FECrsGraph ff_graph(Copy, face_map, 11);
  for (int j = 0; j < cell_map.NumMyElements(); ++j) {
    // Get the cell face indices; we need both process-local and global indices.
    mesh->cell_to_faces((unsigned int) j, (unsigned int*) l_indices, (unsigned int*) l_indices+6);
    for (int i = 0; i < 6; ++i) g_indices[i] = face_map_use.GID(l_indices[i]);
    cf_graph.InsertMyIndices(j, 6, l_indices);
    ff_graph.InsertGlobalIndices(6, g_indices, 6, g_indices);
  }
  cf_graph.FillComplete(face_map, cell_map);
  ff_graph.GlobalAssemble();

  // Create the diagonal cell-cell matrix.
  Dcc_ = new Epetra_Vector(cell_map);

  // Create the cell-face matrix; structure only, no values yet.
  Dcf_ = new Epetra_CrsMatrix(Copy, cf_graph);

  // Create the face-face matrix; structure only, no values yet.
  Dff_ = new Epetra_FECrsMatrix(Copy, ff_graph);
  Dff_->GlobalAssemble();  // need FillComplete() to be true for later copy

  // The face-face Schur complement matrix will be created later if needed.
  Sff_ = 0;
}


DiffusionMatrix::~DiffusionMatrix()
{
  if (Sff_) delete Sff_;
  delete Dff_, Dcf_, Dcc_;
}


const Epetra_FECrsMatrix& DiffusionMatrix::Sff()
{
  // We delay creating the array until requested.
  if (!Sff_) Sff_ = new Epetra_FECrsMatrix(*Dff_);
  return *Sff_;
}


void DiffusionMatrix::Compute(const std::vector<double> &K)
{
  Compute<double>(K);
}


void DiffusionMatrix::Compute(const std::vector<Epetra_SerialSymDenseMatrix> &K)
{
  Compute<Epetra_SerialSymDenseMatrix>(K);
}


template <typename T>
void DiffusionMatrix::Compute(const std::vector<T> &K)
{
  Epetra_Map cell_map = CellMap(false); // owned cells
  Epetra_Map face_map_use = FaceMap(true);  // all used faces

  int l_indices[6];
  Epetra_IntSerialDenseVector g_indices(6);
  Epetra_SerialDenseMatrix minv(6,6);

  (*Dff_).PutScalar(0.0);
  for (int j = 0; j < cell_map.NumMyElements(); ++j) {

    // Compute the inverse of the cell face mass matrix.
    double x[24];
    mesh_->cell_to_coordinates((unsigned int) j, x, x+24);
    MimeticHexLocal mhex(x);
    mhex.mass_matrix(minv, K[j], true);

    // Get the cell face indices; we need both process-local and global indices.
    mesh_->cell_to_faces((unsigned int) j, (unsigned int*) l_indices, (unsigned int*) l_indices+6);
    for (int k = 0; k < 6; ++k) g_indices[k] = face_map_use.GID(l_indices[k]);

    double w[6];
    double matsum = 0.0;
    for (int k = 0; k < 6; ++k) {
      double colsum = 0.0;
      for (int i = 0; i < 6; ++i) colsum += minv(i,k);
      w[k] = -colsum;
      matsum += colsum;
    }

    (*Dcc_)[j] = matsum;
    (*Dcf_).ReplaceMyValues(j, 6, w, l_indices);
    (*Dff_).SumIntoGlobalValues(g_indices, minv);
  }
  ApplyDirichletProjection(*Dff_);
  (*Dff_).GlobalAssemble();
}


void DiffusionMatrix::ComputeFaceSchur()
{
  // Copy matrix values from Dff.
  if (!Sff_) {
    Sff_ = new Epetra_FECrsMatrix(*Dff_);
  } else {
    // Copy matrix values from Dff into Sff; we don't just
    // do Sff = Dff because it does other stuff we don't want.
    int *offsets, *indices;
    double *source, *target;
    (*Dff_).ExtractCrsDataPointers(offsets, indices, source);
    (*Sff_).ExtractCrsDataPointers(offsets, indices, target);
    for (int j = 0; j < (*Dff_).NumMyNonzeros(); ++j) target[j] = source[j];
  }

  double *w;
  int n, *l_indices;
  for (int j = 0; j < (*Dcc_).Map().NumMyElements(); ++j) {
    (*Dcf_).ExtractMyRowView(j, n, w, l_indices);
    Epetra_SerialDenseMatrix update(n,n);
    for (int k = 0; k < n; ++k) {
      for (int i = 0; i < n; ++i) {
        update(i,k) = -w[i]*(w[k]/(*Dcc_)[j]);
      }
    }
    Epetra_IntSerialDenseVector g_indices(n);
    for (int k = 0; k < n; ++k) g_indices[k] = (*Dcf_).ColMap().GID(l_indices[k]);
    (*Sff_).SumIntoGlobalValues(g_indices, update);
  }
  ApplyDirichletProjection(*Sff_);
  (*Sff_).GlobalAssemble();
}

// For each Dirichlet-type face, we need to zero-out the corresponding row and
// column of the face-face CRS matrix and put a 1 on the diagonal. Dealing with
// the row is easy, but zeroing out the column efficiently in parallel is
// tricky and made very difficult by the use of Epetra_FECrsMatrix, which hides
// the critical overlap info needed to do it efficiently.  The approach taken
// here overwrites the values before the parallel assembly of the matrix, and it
// relies on the particular character of this matrix and is not generally valid.
//
// A Dirichlet-type face, being a boundary face, belongs to a single cell, and
// the process that owns the face also owns the cell containing it (a reasonable
// constraint on the distributed mesh).  The non-zero elements of the matrix
// are precisely those that correspond to pairs of faces belonging to a common
// cell.  Thus the non-zeros in the row/column of a Dirichlet-type face are
// precisely those corresponding to the faces of the single cell that contains
// the Dirchlet-type face.  If a process owns a Dirichlet-type face, the values
// in its row/column were computed on the same process when the cell containing
// the face (also owned by the process) was visited.  So what we need to do here
// is overwrite those values with 0.
//
// This function must be called before calling Assemble() on the matrix.

void DiffusionMatrix::ApplyDirichletProjection(Epetra_FECrsMatrix &Mff) const
{
  const Epetra_Map &face_map = FaceMap(false);  // owned faces
  const Epetra_CrsGraph &graph = Mff.Graph();

  int n, *indices;
  const double ZERO = 0.0, ONE  = 1.0;

  for (int j = 0; j < dir_faces_.size(); ++j) {
    int lrow = dir_faces_[j];
    if (!graph.RowMap().MyLID(lrow)) continue; // not my Dirichlet face
    graph.ExtractMyRowView(lrow, n, indices);
    int grow = graph.RowMap().GID(lrow);
    for (int i = 0; i < n; ++i) {
      int gcol = graph.ColMap().GID(indices[i]);
      if (gcol == grow)
        Mff.ReplaceGlobalValues(grow, 1, &ONE, &gcol); // diagonal
      else {
        Mff.ReplaceGlobalValues(grow, 1, &ZERO, &gcol); // row element
        Mff.ReplaceGlobalValues(gcol, 1, &ZERO, &grow); // col element
      }
    }
  }
}

// Sets vector components corresponding Dirichlet-type faces to 0.
// The vector may be based on either the owned or used face map.
// Note that dir_faces_ is the list of all used Dirichlet faces.

void DiffusionMatrix::ApplyDirichletProjection(Epetra_MultiVector &xf) const
{
  for (int j = 0; j < dir_faces_.size(); ++j) {
    if (xf.Map().MyLID(dir_faces_[j]))
      for (int i = 0; i < xf.NumVectors(); ++i) xf[i][dir_faces_[j]] = 0.0;
  }
}


void DiffusionMatrix::Print(std::ostream &os) const
{
  os << *Dcc_ << std::endl;
  os << *Dcf_ << std::endl;
  os << *Dff_ << std::endl;
  if (Sff_) os << *Sff_ << std::endl;
  for (int j = 0; j < dir_faces_.size(); ++j) os << dir_faces_[j] << " ";
  if (dir_faces_.size() > 0) os << std::endl;
}

