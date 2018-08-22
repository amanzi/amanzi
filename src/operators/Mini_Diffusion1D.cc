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

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialize 1D uniform mesh with given end-point areas.
****************************************************************** */
void Mini_Diffusion1D::Init(
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

  diag_.Reshape(ncells + 1);
  up_.Reshape(ncells - 1);
}


/* ******************************************************************
* Create a tri-diagonal matrix
****************************************************************** */
void Mini_Diffusion1D::UpdateMatrices()
{
  int ncells = mesh_->NumRows() - 1; 
  double al, ar, hl, hr, area, x0, x1;

  const auto& mesh = *mesh_;

  x0 = mesh(0);
  x1 = mesh(ncells);
  
  hl = (*K_)(0) / (mesh(1) - mesh(0));
  al = 2 * hl * area_min_;

  for (int i = 0; i < ncells - 1; ++i) {
    double x = mesh(i + 1);
    area = (area_max_ * (mesh(i) - x0) + area_min_ * (x1 - x)) / (x1 - x0);
    area = std::pow(area, igeo_);

    hr = (*K_)(i + 1) / (mesh(i + 2) - x);
    ar = 2 * hl * hr / (hl + hr) * area;

    diag_(i) = al + ar;
    up_(i) = -ar;

    al = ar;
  }

  diag_(ncells) = al;
}

}  // namespace Operators
}  // namespace Amanzi



