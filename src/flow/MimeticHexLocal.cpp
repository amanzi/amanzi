#include <iostream>

#include "cell_geometry.hh"
#include "cell_topology.hh"
using namespace cell_topology;

#include "MimeticHexLocal.hpp"

#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_SerialSpdDenseSolver.h"

void MimeticHexLocal::update(const Epetra_SerialDenseMatrix &x)
{
  cell_geometry::compute_hex_volumes(x, hvol, cwgt);
  double sum_cvol = 0.0;
  for (int i = 0; i < 8; ++i) sum_cvol += cwgt[i];
  for (int i = 0; i < 8; ++i) cwgt[i] = cwgt[i] / sum_cvol;

  //compute_hex_face_normals();
  face_normal.Shape(3,6);
  cell_geometry::compute_hex_face_normals(x, face_normal);
}


void MimeticHexLocal::update(double x[][3])
{
  Epetra_SerialDenseMatrix X(View, (double*)x, 3, 3, 8);
  update(X);
}


void MimeticHexLocal::update(double *x)
{
  Epetra_SerialDenseMatrix X(View, x, 3, 3, 8);
  update(X);
}


void MimeticHexLocal::mass_matrix(Epetra_SerialDenseMatrix &matrix, double K, bool invert) const
{
  Epetra_SerialDenseMatrix Nc(3,3); // face normals at a corner
  Epetra_SerialDenseMatrix Mc(3,3); // corner mass matrix
  Epetra_SerialSymDenseMatrix Mc_sym_view(View, Mc.A(), Mc.LDA(), Mc.M());  // symmetric view of Mc

  // Set the mass matrix to zero.
  for (int j = 0; j < 6; j++)
    for (int i = 0; i < 6; ++i)
      matrix(i,j) = 0.0;

  // Accumulate the mass matrix contribution from each of the eight corners.
  for (int c = 0; c < 8; ++c) {

    // Gather the three face normals adjacent to the corner.
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
	//Nc(k,j) = face_normal[CellTopology::HexCornerFace[c][j]][k];
	Nc(k,j) = face_normal[HexCornerFace[c][j]][k];

    // Mc <-- Nc^T Nc
    Mc.Multiply('T', 'N', 1.0, Nc, Nc, 0.0);

    // Mc <-- Mc^(-1) using the symmetric Choleski method
    Epetra_SerialSpdDenseSolver solver;
    solver.SetMatrix(Mc_sym_view);
    solver.Invert(); // gives full inverse in Mc, not just in triangle

    // Scatter the corner mass matrix into the full cell mass matrix.
    double s = hvol * cwgt[c] / K;
    for (int j = 0; j < 3; ++j) {   // loop over corner face cols
      //int jj = CellTopology::HexCornerFace[c][j]; // the corresponding cell face index
      int jj = HexCornerFace[c][j]; // the corresponding cell face index
      for (int i = 0; i <3 ; ++i) {   // loop over corner face rows
	//int ii = CellTopology::HexCornerFace[c][i]; // the corresponding cell face index
	int ii = HexCornerFace[c][i]; // the corresponding cell face index
        matrix(ii,jj) += s * Mc(i,j);
      }
    }
  }

  if (invert) {
    // Create a "symmetric" view of the matrix (which is symmetric).
    Epetra_SerialSymDenseMatrix sym_view(View, matrix.A(), matrix.LDA(), matrix.M());
    // Create an SPD "solver" object for the matrix.
    Epetra_SerialSpdDenseSolver solver;
    solver.SetMatrix(sym_view);
    // Now invert the matrix in-place using the Choleski method.
    solver.Invert(); // gives full inverse in matrix, not just in triangle.
  }
}


void MimeticHexLocal::mass_matrix(Epetra_SerialDenseMatrix &matrix, const Epetra_SerialSymDenseMatrix &K, bool invert) const
{
  Epetra_SerialDenseMatrix Nc(3,3); // face normals at a corner
  Epetra_SerialDenseMatrix Mc(3,3); // corner mass matrix
  Epetra_SerialDenseMatrix Tc(3,3); // temporary corner matrix
  Epetra_SerialSymDenseMatrix Mc_sym_view(View, Mc.A(), Mc.LDA(), Mc.M());  // symmetric view of Mc

  // Set the mass matrix to zero.
  for (int j = 0; j < 6; j++)
    for (int i = 0; i < 6; ++i)
      matrix(i,j) = 0.0;

  // Accumulate the mass matrix contribution from each of the eight corners.
  for (int c = 0; c < 8; ++c) {

    // Gather the three face normals adjacent to the corner.
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
	//Nc(k,j) = face_normal[CellTopology::HexCornerFace[c][j]][k];
	Nc(k,j) = face_normal[HexCornerFace[c][j]][k];

    // Mc <-- Nc^T K Nc
    Tc.Multiply('L', 1.0, K, Nc, 0.0);
    Mc.Multiply('T', 'N', 1.0, Nc, Tc, 0.0);

    // Mc <-- Mc^(-1) using the symmetric Choleski method
    Epetra_SerialSpdDenseSolver solver;
    solver.SetMatrix(Mc_sym_view);
    solver.Invert(); // gives full inverse in Mc, not just in triangle

    // Scatter the corner mass matrix into the full cell mass matrix.
    double s = hvol * cwgt[c];
    for (int j = 0; j < 3; ++j) {   // loop over corner face cols
      //int jj = CellTopology::HexCornerFace[c][j]; // the corresponding cell face index
      int jj = HexCornerFace[c][j]; // the corresponding cell face index
      for (int i = 0; i <3 ; ++i) {   // loop over corner face rows
	//int ii = CellTopology::HexCornerFace[c][i]; // the corresponding cell face index
	int ii = HexCornerFace[c][i]; // the corresponding cell face index
        matrix(ii,jj) += s * Mc(i,j);
      }
    }
  }

  if (invert) {
    // Create a "symmetric" view of the matrix (which is symmetric).
    Epetra_SerialSymDenseMatrix sym_view(View, matrix.A(), matrix.LDA(), matrix.M());
    // Create an SPD "solver" object for the matrix.
    Epetra_SerialSpdDenseSolver solver;
    solver.SetMatrix(sym_view);
    // Now invert the matrix in-place using the Choleski method.
    solver.Invert(); // gives full inverse in matrix, not just in triangle.
  }
}


void MimeticHexLocal::diff_op(double coef,
    const double &pcell, const double pface[],
    double &rcell, double rface[]) const
{
  Epetra_SerialDenseVector aux1(6);
  Epetra_SerialDenseMatrix Minv(6,6);

  // Inverse of the mass matrix.
  mass_matrix(Minv, coef, true);

  for (int i = 0; i < 6; ++i)
    aux1(i) = pface[i] - pcell;

  Epetra_SerialDenseVector aux2(View, rface, 6);

  Minv.Multiply(false, aux1, aux2);

  rcell = 0.0;
  for (int i = 0; i < 6; ++i)
    rcell -= rface[i];
}


void MimeticHexLocal::diff_op(const Epetra_SerialSymDenseMatrix &coef,
    const double &pcell, const double pface[],
    double &rcell, double rface[]) const
{
  Epetra_SerialDenseVector aux1(6);
  Epetra_SerialDenseMatrix Minv(6,6);

  // Inverse of the mass matrix.
  mass_matrix(Minv, coef, true);

  for (int i = 0; i < 6; ++i)
    aux1(i) = pface[i] - pcell;

  Epetra_SerialDenseVector aux2(View, rface, 6);

  Minv.Multiply(false, aux1, aux2);

  rcell = 0.0;
  for (int i = 0; i < 6; ++i)
    rcell -= rface[i];
}


void MimeticHexLocal::diff_op(double coef,
    const double &pcell, const Epetra_SerialDenseVector &pface,
    double &rcell, Epetra_SerialDenseVector &rface) const
{
  Epetra_SerialDenseVector aux(6);
  Epetra_SerialDenseMatrix Minv(6,6);

  // Inverse of the mass matrix.
  mass_matrix(Minv, coef, true);

  for (int i = 0; i < 6; ++i)
    aux(i) = pface[i] - pcell;

  Minv.Multiply(false, aux, rface);

  rcell = 0.0;
  for (int i = 0; i < 6; ++i)
    rcell -= rface(i);
}


void MimeticHexLocal::diff_op(const Epetra_SerialSymDenseMatrix &coef,
    const double &pcell, const Epetra_SerialDenseVector &pface,
    double &rcell, Epetra_SerialDenseVector &rface) const
{
  Epetra_SerialDenseVector aux(6);
  Epetra_SerialDenseMatrix Minv(6,6);

  // Inverse of the mass matrix.
  mass_matrix(Minv, coef, true);

  for (int i = 0; i < 6; ++i)
    aux(i) = pface[i] - pcell;

  Minv.Multiply(false, aux, rface);

  rcell = 0.0;
  for (int i = 0; i < 6; ++i)
    rcell -= rface(i);
}

//void MimeticHexLocal::Flux(double &coef, const double &pcell, const Epetra_SerialDenseVector &pface, Epetra_SerialDenseVector &fface) const
//{
//  mass_matrix(Minv, coef, true);
//  for (int i = 0; i < 6; ++i) {
//    aux[i] = pcell - pface[i];
//    Minv.Multiply(false, aux, fface);
//  }
//}


void MimeticHexLocal::GravityFlux(const double g[], double gflux[]) const
{
  // Expect g[3] and gflux[6]
  for (int i = 0; i < 6; ++i) { // loop over faces
    double s = 0.0;
    for (int k = 0; k < 3; ++k)
      s += g[k] * face_normal[i][k];
    gflux[i] = s;
  }
}


void MimeticHexLocal::CellFluxVector(double Fface[], double Fcell[]) const
{
  Epetra_SerialSymDenseMatrix a;
  Epetra_SerialDenseVector b;
  a.Shape(3);
  b.Size(3);
  for (int k = 0; k < 6; ++k) {
    double w = 1.0 / cell_geometry::vector_length(face_normal[k],3);
    //double w = 1.0;
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i)
        a(i,j) += w * face_normal(i,k) * face_normal(j,k);
      b(j) += w * Fface[k] * face_normal(j,k);
    }
  }
  Epetra_SerialSpdDenseSolver solver;
  solver.SetMatrix(a);
  Epetra_SerialDenseVector x;
  x.Size(3);
  solver.SetVectors(x, b);
  solver.Factor();
  solver.Solve();
  for (int k = 0; k < 3; ++k) Fcell[k] = x[k];
}

