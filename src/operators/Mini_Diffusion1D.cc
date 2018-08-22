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

#include <Mini_Diffusion1D.hh>

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Simple constructor
****************************************************************** */
Mini_Diffusion1D::Mini_Diffusion1D(std::shared_ptr<const WhetStone::DenseVector> mesh)
  : mesh_(mesh)
{
  int ncells = mesh_->NumRows();
  diag_.Reshape(ncells - 1);
}


/* ******************************************************************
* Create a tri-diagonal matrix
****************************************************************** */
void Mini_Diffusion1D::UpdateMatrices()
{
  int ncells = mesh_->NumRows();
  for (int i = 0; i < ncells; ++i) {
    diag_(i) = (*K_)(i) + (*K_)(i + 1);
  }
}

}  // namespace Operators
}  // namespace Amanzi



