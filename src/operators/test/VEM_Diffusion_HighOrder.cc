/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Virtual element high-order conformal discretization of elliptic operator
  via a static condensation.
*/

#include <cmath>
#include <iterator>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "VEM_Diffusion_HighOrder.hh"
#include "Tensor.hh"
#include "VEM_RaviartThomasSerendipity.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

RegisteredFactory<VEM_Diffusion_HighOrder> VEM_Diffusion_HighOrder::factory_("vem diffusion high order");

/* ******************************************************************
* Stiffness matrix via static condensation
****************************************************************** */
int VEM_Diffusion_HighOrder::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix M;
  VEM_RaviartThomasSerendipity vem(plist_, mesh_);

  int ok = vem.MassMatrix(c, K, M);
  if (ok) return ok;

  A = M;
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi



