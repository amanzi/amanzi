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
  plist_.set<bool>("save face matrices", true);
  VEM_RaviartThomasSerendipity vem(plist_, mesh_);

  int ok = vem.MassMatrix(c, K, M);
  if (ok) return ok;

  // static condensation
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const auto& MGf = vem.MGf(); 
  AMANZI_ASSERT(MGf.size() == nfaces);

  int nrows = M.NumRows();
  DenseMatrix C(nrows, nrows), B(nrows, 1), BT(1, nrows);
  C.PutScalar(0.0);
  B.PutScalar(0.0);
  BT.PutScalar(0.0);

  // -- populate interface and divergence matrices
  int pos(0);
  for (int n = 0; n < nfaces; ++n) {
    double area = mesh_->face_area(faces[n]);

    int mrows = MGf[n].NumRows(); 
    C.InsertSubMatrix(MGf[n], 0, mrows, 0, mrows, pos, pos);
    B(pos, 0) = BT(0, pos) = area;

    pos += mrows;
  }
  
  // -- elliminate lagrange multipliers, leave central pressure
  M.InverseSPD();
  auto MC = M * C;
  auto MB = M * B;
  
  auto CMC = C * MC;
  auto BMC = BT * MC;
  auto CMB = C * MB;
  auto BMB = BT * MB;

  A.Reshape(nrows + 1, nrows + 1);
  A.InsertSubMatrix(CMC, 0, nrows, 0, nrows, 0, 0);
  A.InsertSubMatrix(CMB, 0, nrows, 0, 1, 0, nrows);
  A.InsertSubMatrix(BMC, 0, 1, 0, nrows, nrows, 0);
  A(nrows, nrows) = BMB(0, 0);

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi



