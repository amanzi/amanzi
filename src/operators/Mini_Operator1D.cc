/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  This provides a base class for other mini classes that implement 
  mathematical models for special physics.
*/

#include <dbc.hh>

#include <Mini_Operator1D.hh>

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Make a hard copy
****************************************************************** */
Mini_Operator1D::Mini_Operator1D(const Mini_Operator1D& other)
  : mesh_(other.mesh_),
    diag_(other.diag_),
    up_(other.up_),
    down_(other.down_),
    rhs_(other.rhs_)
{}


/* ******************************************************************
* Initialize 1D uniform mesh with given end-point areas and allocate
* matrix for FV-type discretizations.
****************************************************************** */
void Mini_Operator1D::Init(std::shared_ptr<const WhetStone::DenseVector> mesh)
{
  mesh_ = mesh;

  int ncells = mesh_->NumRows() - 1;
  AMANZI_ASSERT(ncells > 0);

  diag_.Reshape(ncells);
  up_.Reshape(ncells);
  down_.Reshape(ncells);

  rhs_.Reshape(ncells);
  rhs_.PutScalar(0.0);
}


/* ******************************************************************
* Update matrix and right-hand side: linear model
****************************************************************** */
void Mini_Operator1D::AddAccumulationTerm(double s0, double s1, double dt,
                                          WhetStone::DenseVector& sol) 
{
  int ncells = diag_.NumRows();
  for (int i = 0; i < ncells; ++i) {
    double h = (*mesh_)(i + 1) - (*mesh_)(i);
    diag_(i) += s1 * h / dt;
    rhs_(i) += s0 * sol(i) * h / dt;
  }
}

void Mini_Operator1D::AddAccumulationTerm(
    const WhetStone::DenseVector& s0, const WhetStone::DenseVector& s1,
    double dt, WhetStone::DenseVector& sol) 
{
  int ncells = diag_.NumRows();
  for (int i = 0; i < ncells; ++i) {
    double h = mesh_cell_volume(i);
    diag_(i) += s1(i) * h / dt;
    rhs_(i) += s0(i) * sol(i) * h / dt;
  }
}

void Mini_Operator1D::AddAccumulationTerm(const WhetStone::DenseVector& s1, bool add_volume)
{
  int ncells = diag_.NumRows();
  if (add_volume) {
    for (int i = 0; i < ncells; ++i) {
      diag_(i) += s1(i) * mesh_cell_volume(i);
    }
  } else {
    diag_ += s1;
  }
}


/* ******************************************************************
* Matrix-vector product
****************************************************************** */
void Mini_Operator1D::Apply(const WhetStone::DenseVector& v,
                            WhetStone::DenseVector& av)
{
  int ncells = diag_.NumRows();
  for (int i = 0; i < ncells; ++i) {
    av(i) = diag_(i) * v(i);
  }
  for (int i = 0; i < ncells - 1; ++i) {
    av(i) += up_(i) * v(i + 1);
  }
  for (int i = 1; i < ncells; ++i) {
    av(i) += down_(i) * v(i - 1);
  }
}


/* ******************************************************************
* Solve linear system using direct method
****************************************************************** */
void Mini_Operator1D::ApplyInverse(const WhetStone::DenseVector& rhs,
                                   WhetStone::DenseVector& sol)
{
  int nrhs(1), info;
  int n = diag_.NumRows();

  double* dl = down_.Values(); 
  double* dr = up_.Values(); 
  dl++;

  sol = rhs;
  WhetStone::DGTSV_F77(&n, &nrhs, dl, diag_.Values(), dr,
                       sol.Values(), &n, &info);
}


/* ******************************************************************
* Matrix modifications
****************************************************************** */
void Mini_Operator1D::GetMatrixRow(int i, double* al, double* ad, double* ar) const
{
  *al = down_(i);
  *ad = diag_(i);
  *ar = up_(i);
}

void Mini_Operator1D::SetMatrixRow(int i, double al, double ad, double ar)
{
  down_(i) = al;
  diag_(i) = ad;
  up_(i) = ar;
}


/* ******************************************************************
* Visualize matrix
****************************************************************** */
void Mini_Operator1D::Print(int n, const char* format)
{
  double zero(0.0);
  int m = (n > 0) ? n : diag_.NumRows();

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < i - 1; ++j) printf(format, zero);
    if (i > 0) printf(format, down_(i));
    printf(format, diag_(i));
    if (i < diag_.NumRows() - 1) printf(format, up_(i));
    for (int j = i + 1; j < m; ++j) printf(format, zero);
    printf("\n");
  }
}

}  // namespace Operators
}  // namespace Amanzi



