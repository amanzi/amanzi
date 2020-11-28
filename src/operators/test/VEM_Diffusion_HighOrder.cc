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
  const auto& faces = mesh_->cell_get_faces(c);
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
    B(pos, 0) = BT(0, pos) = -area;  // trick to make off-diagonal blocks negative

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

  return 0;
}


/* ******************************************************************
* Post-Processing
****************************************************************** */
void VEM_Diffusion_HighOrder::UpdateFlux(
    int c, const Tensor& K, double** lambda, const double* p, double** flux)
{
  // calculate local matrices (costly, good only for test)
  DenseMatrix M;
  plist_.set<bool>("save face matrices", true);
  VEM_RaviartThomasSerendipity vem(plist_, mesh_);

  vem.MassMatrix(c, K, M);
  int nrows = M.NumRows(); 

  // calculate flux, face-by-face
  const auto& faces = mesh_->cell_get_faces(c);
  const auto& dirs = mesh_->cell_get_face_dirs(c);
  int nfaces = faces.size();

  const auto& MGf = vem.MGf(); 
  int ndf = MGf[0].NumRows();

  DenseVector dp(nrows), uf(nrows), ploc(ndf), tmp(ndf);

  int pos(0);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    double area = mesh_->face_area(f);

    for (int i = 0; i < ndf; ++i) ploc(i) = lambda[i][f];
    MGf[n].Multiply(ploc, tmp, false);
    for (int i = 0; i < ndf; ++i) dp(pos + i) = -tmp(i);
    dp(pos) += area * p[c];

    pos += ndf;
  }

  M.InverseSPD();
  M.Multiply(dp, uf, false);

  // copy flux to the global vector and scale it by area
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    double area = mesh_->face_area(f);

    for (int i = 0; i < ndf; ++i)
      flux[i][f] = uf(ndf * n + i) * dirs[n] * area;
  }
}

}  // namespace WhetStone
}  // namespace Amanzi



