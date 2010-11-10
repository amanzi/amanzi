#include "DiffusionMatrix.hpp"

#include "Epetra_FECrsGraph.h"
#include "MimeticHexLocal.hpp"

#include <iostream>

DiffusionMatrix::DiffusionMatrix(const Teuchos::RCP<Mesh_maps_base> &mesh, const std::vector<int> &dir_faces)
{
  mesh_ = mesh;
  dir_faces_ = dir_faces;
  Sff_ = 0;

  Epetra_Map cell_map = mesh->cell_map(false);
  Epetra_Map face_map = mesh->face_map(false);
  Epetra_Map face_map_ovl = mesh->face_map(true);

  Dcc_ = new Epetra_Vector(cell_map);

  //
  // CREATE THE CELL-FACE MATRIX

  unsigned int cface[6];
  int indices[6];

  Epetra_CrsGraph cf_graph(Copy, cell_map, face_map_ovl, 6, true);
  for (unsigned int j = cell_map.MinLID(); j <= cell_map.MaxLID(); ++j) {
    mesh->cell_to_faces(j, cface, cface+6); // process-local face indices
    for (int i = 0; i < 6; ++i) indices[i] = cface[i]; // just convert to signed int
    cf_graph.InsertMyIndices(j, 6, indices);
  }
  cf_graph.FillComplete(face_map_ovl,cell_map);
  Dcf_ = new Epetra_CrsMatrix(Copy, cf_graph);

  //
  // CREATE THE FACE-FACE MATRIX GRAPH

  Epetra_FECrsGraph graph(Copy, face_map, 11);
  for (unsigned int j = cell_map.MinLID(); j <= cell_map.MaxLID(); ++j) {
    mesh->cell_to_faces(j, cface, cface+6); // process-local face indices, possibly ghosts
    for (int i = 0; i < 6; ++i) indices[i] = face_map_ovl.GID(cface[i]); // global face indices
    graph.InsertGlobalIndices(6, indices, 6, indices);
  }
  graph.GlobalAssemble();
  Dff_ = new Epetra_FECrsMatrix(Copy, graph);
  Dff_->GlobalAssemble();  // need FillComplete() to be true for later copy
}


DiffusionMatrix::~DiffusionMatrix()
{
  delete Dcc_, Dcf_, Dff_;
  if (Sff_) delete Sff_;
}


void DiffusionMatrix::Compute(const std::vector<double> &K)
{
  // Local storage for the 8 vertex coordinates of a hexahedral cell.
  double x[8][3];
  double *xBegin = &x[0][0];  // begin iterator
  double *xEnd = xBegin+24;   // end iterator

  Epetra_SerialDenseMatrix minv(6,6);
  const bool invert = true;

  unsigned int cface[6];
  int l_indices[6];
  Epetra_IntSerialDenseVector g_indices(6);

  Epetra_Map cell_map = mesh_->cell_map(false);
  Epetra_Map face_map = mesh_->face_map(true);

  (*Dff_).PutScalar(0.0);
  for (int icell = 0; icell <= cell_map.MaxLID(); ++icell) {

    // Compute the inverse of the cell face mass matrix.
    mesh_->cell_to_coordinates((unsigned int) icell, xBegin, xEnd);
    MimeticHexLocal mhex(x);
    mhex.mass_matrix(minv, K[icell], invert);
\
    // Get the cell face indices; we need both local and global indices.
    // mesh_->cell_to_faces((unsigned int) icell, l_indices, l_indices+6);
    mesh_->cell_to_faces((unsigned int) icell, cface, cface+6);
    for (int k = 0; k < 6; ++k) {
      l_indices[k] = (int) cface[k]; // need ints, not unsigned ints
      g_indices[k] = face_map.GID(l_indices[k]);
    }

    double w[6];
    double matsum = 0.0;
    for (int j = 0; j < 6; ++j) {
      double colsum = 0.0;
      for (int i = 0; i < 6; ++i) {
        colsum += minv(i,j);
      }
      w[j] = -colsum;
      matsum += colsum;
    }

    (*Dcc_)[icell] = matsum;
    (*Dcf_).ReplaceMyValues(icell, 6, w, l_indices);
    (*Dff_).SumIntoGlobalValues(g_indices, minv);
  }
  (*Dcf_).FillComplete(face_map,cell_map);
  (*Dff_).GlobalAssemble();

  ApplyDirichletProjection(*Dff_);
}


void DiffusionMatrix::Compute(const std::vector<Epetra_SerialSymDenseMatrix> &K)
{
  // Local storage for the 8 vertex coordinates of a hexahedral cell.
  double x[8][3];
  double *xBegin = &x[0][0];  // begin iterator
  double *xEnd = xBegin+24;   // end iterator

  Epetra_SerialDenseMatrix minv(6,6);
  const bool invert = true;

  unsigned int cface[6];
  int l_indices[6];
  Epetra_IntSerialDenseVector g_indices(6);

  Epetra_Map cell_map = mesh_->cell_map(false);
  Epetra_Map face_map = mesh_->face_map(true);

  (*Dff_).PutScalar(0.0);
  for (int icell = 0; icell <= cell_map.MaxLID(); ++icell) {

    // Compute the inverse of the cell face mass matrix.
    mesh_->cell_to_coordinates((unsigned int) icell, xBegin, xEnd);
    MimeticHexLocal mhex(x);
    mhex.mass_matrix(minv, K[icell], invert);
\
    // Get the cell face indices; we need both local and global indices.
    // mesh_->cell_to_faces((unsigned int) icell, l_indices, l_indices+6);
    mesh_->cell_to_faces((unsigned int) icell, cface, cface+6);
    for (int k = 0; k < 6; ++k) {
      l_indices[k] = (int) cface[k]; // need ints, not unsigned ints
      g_indices[k] = face_map.GID(l_indices[k]);
    }

    double w[6];
    double matsum = 0.0;
    for (int j = 0; j < 6; ++j) {
      double colsum = 0.0;
      for (int i = 0; i < 6; ++i) {
        colsum += minv(i,j);
      }
      w[j] = -colsum;
      matsum += colsum;
    }

    (*Dcc_)[icell] = matsum;
    (*Dcf_).ReplaceMyValues(icell, 6, w, l_indices);
    (*Dff_).SumIntoGlobalValues(g_indices, minv);
  }
  (*Dcf_).FillComplete(face_map,cell_map);
  (*Dff_).GlobalAssemble();

  ApplyDirichletProjection(*Dff_);
}


const Epetra_FECrsMatrix& DiffusionMatrix::Sff()
{
  // We delay creating the array until requested.
  if (!Sff_) Sff_ = new Epetra_FECrsMatrix(*Dff_);
  return *Sff_;
}


void DiffusionMatrix::ComputeFaceSchur()
{
  // Copy matrix values from Dff.
  if (!Sff_) {
    Sff_ = new Epetra_FECrsMatrix(*Dff_);
  } else {
    // Copy matrix values from Dff into Sff
    int *offsets, *indices;
    double *source, *target;
    (*Dff_).ExtractCrsDataPointers(offsets, indices, source);
    (*Sff_).ExtractCrsDataPointers(offsets, indices, target);
    for (int j = 0; j < (*Dff_).NumMyNonzeros(); ++j)
      target[j] = source[j];
  }

  int n, *indices;
  double *w;
  for (int icell = 0; icell <= (*Dcc_).Map().MaxLID(); ++icell) {
    (*Dcf_).ExtractMyRowView(icell, n, w, indices);
    Epetra_SerialDenseMatrix update(n,n);
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        update(i,j) = -w[i]*(w[j]/(*Dcc_)[icell]);
      }
    }
    for (int i = 0; i < n; ++i) indices[i] = (*Dcf_).ColMap().GID(indices[i]);
    Epetra_IntSerialDenseVector g_indices(View, indices, n);
    (*Sff_).SumIntoGlobalValues(g_indices, update);
  }
  (*Sff_).GlobalAssemble();

  ApplyDirichletProjection(*Sff_);
}


void DiffusionMatrix::ApplyDirichletProjection(Epetra_CrsMatrix &Mff)
{
  int n;
  int *indices;
  double *values;
  const double ZERO = 0.0, ONE  = 1.0;

  for (int j = 0; j < dir_faces_.size(); ++j) {
    int lrow = dir_faces_[j]; // local row index
    Mff.ExtractMyRowView(lrow, n, values, indices);
    for (int i = 0; i < n; ++i) {
      values[i] = 0.0;
      int lcol = indices[i]; // local column index
      if (Mff.RowMap().MyLID(lcol)) Mff.ReplaceMyValues(lcol, 1, &ZERO, &lrow);
    }
    Mff.ReplaceMyValues(lrow, 1, &ONE, &lrow); // put a 1 on the diagonal
  }
}


void DiffusionMatrix::ApplyDirichletProjection(Epetra_MultiVector &xf)
{
  for (int j = 0; j < dir_faces_.size(); ++j)
    for (int i = 0; i < xf.NumVectors(); ++i)
      xf[i][dir_faces_[j]] = 0.0;
}


void DiffusionMatrix::Print()
{
  std::cout << *Dcc_ << std::endl;
  std::cout << *Dcf_ << std::endl;
  std::cout << *Dff_ << std::endl;
  if (Sff_) std::cout << *Sff_ << std::endl;
}

