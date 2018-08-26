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
void Mini_Diffusion1D::UpdateMatrices()
{
  int ncells = mesh_->NumRows() - 1; 
  double al, ar, hl, hr, area, x0, x1, Kc;

  const auto& mesh = *mesh_;

  x0 = mesh(0);
  x1 = mesh(ncells);
  
  Kc = (K_ != NULL) ? (*K_)(0) : Kconst_;
  hl = Kc / (mesh(1) - mesh(0));
  al = 2 * hl * area_min_;

  for (int i = 0; i < ncells - 1; ++i) {
    double x = mesh(i + 1);
    area = (area_max_ * (x - x0) + area_min_ * (x1 - x)) / (x1 - x0);
    area = std::pow(area, igeo_);

    Kc = (K_ != NULL) ? (*K_)(i + 1) : Kconst_;
    hr = Kc / (mesh(i + 2) - x);
    ar = 2 * hl * hr / (hl + hr) * area;

    diag_(i) = al + ar;
    down_(i) = -al;
    up_(i) = -ar;

    al = ar;
    hl = hr;
  }

  Kc = (K_ != NULL) ? (*K_)(ncells - 1) : Kconst_;
  hr = Kc / (mesh(ncells) - mesh(ncells - 1));
  ar = 2 * hr * area_max_;

  diag_(ncells - 1) = al + ar;
  down_(ncells - 1) = -al;
  up_(ncells - 1) = -ar;
}


/* ******************************************************************
* Apply boundary conditions.
****************************************************************** */
void Mini_Diffusion1D::ApplyBCs(double bcl, int type_l, double bcr, int type_r)
{
  if (type_l == Operators::OPERATOR_BC_DIRICHLET) {
    rhs_(0) -= down_(0) * bcl;
  } else if(type_l == Operators::OPERATOR_BC_NEUMANN) {
    diag_(0) += down_(0);
    rhs_(0) += bcl;
  }

  int n = mesh_->NumRows() - 2;
  if (type_r == Operators::OPERATOR_BC_DIRICHLET) {
    rhs_(n) -= up_(n) * bcr;
  } else if(type_r == Operators::OPERATOR_BC_NEUMANN) {
    diag_(n) += up_(n);
    rhs_(n) += bcr;
  }
}

}  // namespace Operators
}  // namespace Amanzi



