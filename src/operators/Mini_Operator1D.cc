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
* Initialize 1D uniform mesh with given end-point areas and allocate
* matrix for FV-type discretizations.
****************************************************************** */
void Mini_Operator1D::Init(
    std::shared_ptr<const WhetStone::DenseVector> mesh,
    const std::string& geometry, double area_min, double area_max)
{
  mesh_ = mesh;
  area_min_ = area_min;
  area_max_ = area_max;

  int ncells = mesh_->NumRows() - 1;
  AMANZI_ASSERT(ncells > 0);

  igeo_ = 1;
  if (geometry == "spherical") igeo_ = 2;

  diag_.Reshape(ncells);
  up_.Reshape(ncells);
  down_.Reshape(ncells);

  rhs_.Reshape(ncells);
  rhs_.PutScalar(0.0);
}


/* ******************************************************************
* Update matrix and right-hand sides.
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

}  // namespace Operators
}  // namespace Amanzi



