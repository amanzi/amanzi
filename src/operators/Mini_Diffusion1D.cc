/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Mini classes implement mathematical models for special physics, such 
  as serial 1D dual porosity models. 
*/

#include <dbc.hh>

#include <Mini_Diffusion1D.hh>
#include <OperatorDefs.hh>

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Create a tri-diagonal matrix. Structure of underlying equation is
*   down_(i) x_{i-1} + diag_(i) x_i + up_(i) x_{i + 1} = b_i
* We use end values of sub-diagonals to impose boundary conditions.
****************************************************************** */
void Mini_Diffusion1D::UpdateMatrices(const PDEType method)
{
  if (method == PDEType::PDE_DIFFUSION_MFD)
    UpdateMatricesMFD_(); 
  else if (method == PDEType::PDE_DIFFUSION_FD)
    UpdateMatricesFD_(); 
}


/* ******************************************************************
* MFD scheme with harmonic mean for K, arithmetic mean for k
****************************************************************** */
void Mini_Diffusion1D::UpdateMatricesMFD_()
{
  int ncells = mesh_->NumRows() - 1; 
  double al, ar, hl, hr, Kc;

  const auto& mesh = *mesh_;

  Kc = (K_ != NULL) ? (*K_)(0) : Kconst_;

  hl = Kc / (mesh(1) - mesh(0));
  al = 2 * hl;
  if (k_ != NULL) al *= (*k_)(0);

  for (int i = 0; i < ncells - 1; ++i) {
    Kc = (K_ != NULL) ? (*K_)(i + 1) : Kconst_;

    hr = Kc / (mesh(i + 2) - mesh(i + 1));
    ar = 2 * hl * hr / (hl + hr);
    if (k_ != NULL) ar *= ((*k_)(i) + (*k_)(i + 1)) / 2;

    diag_(i) = al + ar;
    down_(i) = -al;
    up_(i) = -ar;

    al = ar;
    hl = hr;
  }

  Kc = (K_ != NULL) ? (*K_)(ncells - 1) : Kconst_;

  hr = Kc / (mesh(ncells) - mesh(ncells - 1));
  ar = 2 * hr;
  if (k_ != NULL) ar *= (*k_)(ncells - 1);

  diag_(ncells - 1) = al + ar;
  down_(ncells - 1) = -al;
  up_(ncells - 1) = -ar;
}


/* ******************************************************************
* FV scheme with K=1 and arithmetic for k
****************************************************************** */
void Mini_Diffusion1D::UpdateMatricesFD_()
{
  int ncells = mesh_->NumRows() - 1; 
  double al, ar;

  const auto& mesh = *mesh_;

  al = (*k_)(0) / (mesh(1) - mesh(0));

  for (int i = 0; i < ncells; ++i) {
    ar = (*k_)(i) / (mesh(i + 1) - mesh(i));

    diag_(i) = al + ar;
    down_(i) = -al;
    up_(i) = -ar;

    al = ar;
  }
}


/* ******************************************************************
* Jacobian matrix for operator A(k(p)) p - f(k(p)).
* NOTE: we assume that k_ != NULL and dkdp_ != NULL, i.e. J != A.
****************************************************************** */
void Mini_Diffusion1D::UpdateJacobian(
    const WhetStone::DenseVector& p,
    double bcl, int type_l, double bcr, int type_r)
{
  int ncells = mesh_->NumRows() - 1; 
  double al, ar, bl, br, hl, hr, Kc, tmp0, tmp1;

  const auto& mesh = *mesh_;
  const auto& k = *k_;
  const auto& dkdp = *dkdp_;

  // derivatives of A(k(p))
  Kc = (K_ != NULL) ? (*K_)(0) : Kconst_;
  hl = Kc / mesh_cell_volume(0);
  al = 2 * hl;
  tmp0 = al;
  bl = al * p(0);
  al *= k(0);

  for (int i = 0; i < ncells - 1; ++i) {
    int j = (i == 0) ? 0 : i - 1; 

    Kc = (K_ != NULL) ? (*K_)(i + 1) : Kconst_;
    hr = Kc / (mesh(i + 2) - mesh(i + 1));
    ar = hl * hr / (hl + hr);
    br = ar * (p(i + 1) - p(i));
    ar *= k(i) + k(i + 1);

    diag_(i) = (al + ar) + (bl - br) * dkdp(i);
    down_(i) = -al + bl * dkdp(j);
    up_(i) = -ar - br * dkdp(i + 1);

    al = ar;
    bl = br;
    hl = hr;
  }

  Kc = (K_ != NULL) ? (*K_)(ncells - 1) : Kconst_;
  hr = Kc / mesh_cell_volume(ncells - 1);
  ar = 2 * hr;
  tmp1 = ar;
  br = ar * p(ncells - 1);
  ar *= k(ncells - 1);

  diag_(ncells - 1) = (al + ar) + (bl + br) * dkdp(ncells - 1);
  down_(ncells - 1) = -al + bl * dkdp(ncells - 2);
  up_(ncells - 1) = -ar - br * dkdp(ncells - 1);

  // derivatives of f(k(p))
  if (type_l == Operators::OPERATOR_BC_DIRICHLET) {
    diag_(0) -= tmp0 * dkdp(0) * bcl;
  }

  if (type_r == Operators::OPERATOR_BC_DIRICHLET) {
    diag_(ncells - 1) -= tmp1 * dkdp(ncells - 1) * bcr;
  }
  else if (type_r == Operators::OPERATOR_BC_NEUMANN) {
    diag_(ncells - 1) += up_(ncells - 1);
  }
}


/* ******************************************************************
* Apply boundary conditions.
****************************************************************** */
void Mini_Diffusion1D::ApplyBCs(double bcl, int type_l, double bcr, int type_r)
{
  if (type_l == Operators::OPERATOR_BC_DIRICHLET) {
    rhs_(0) -= down_(0) * bcl;
  } else if (type_l == Operators::OPERATOR_BC_NEUMANN) {
    diag_(0) += down_(0);
    rhs_(0) += bcl;
  }

  int n = mesh_->NumRows() - 2;
  if (type_r == Operators::OPERATOR_BC_DIRICHLET) {
    rhs_(n) -= up_(n) * bcr;
  } else if (type_r == Operators::OPERATOR_BC_NEUMANN) {
    diag_(n) += up_(n);
    rhs_(n) -= bcr;
  }
}

}  // namespace Operators
}  // namespace Amanzi



