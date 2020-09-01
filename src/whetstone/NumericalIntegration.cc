/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Numerical and exact integration over polytopal cells. 
*/

#include <memory>

#include "Teuchos_RCP.hpp"

#include "NumericalIntegration.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Integrate over face f a group of non-normalized monomials of
* the same order k centered at the centroid of cell c.
****************************************************************** */
template <>
void NumericalIntegration<AmanziMesh::Mesh>::IntegrateMonomialsFaceReduction_(
    int c, int f, double factor, int k, Polynomial& integrals) const
{
  int nk = PolynomialSpaceDimension(d_, k - 1);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
  const AmanziGeometry::Point& normal = mesh_->face_normal(f);

  // create a surface mesh
  auto coordsys = std::make_shared<SurfaceCoordinateSystem>(xf, normal);
  Teuchos::RCP<const SurfaceMiniMesh> surf_mesh = Teuchos::rcp(new SurfaceMiniMesh(mesh_, coordsys));
  NumericalIntegration<SurfaceMiniMesh> numi_f(surf_mesh);

  PolynomialIterator it(d_);
  for (it.begin(k); it.MonomialSetOrder() <= k; ++it) {
    int l = it.MonomialSetPosition();

    // using monomial centered at xc, create 2D polynomial centered at xf
    const int* idx = it.multi_index();
    Polynomial poly(d_, idx, 1.0);
    poly.set_origin(xc);
    poly.ChangeCoordinates(xf, *coordsys->tau());  

    integrals(nk + l) += factor * numi_f.IntegratePolynomialCell(f, poly);
  }
}


template <>
void NumericalIntegration<SurfaceMiniMesh>::IntegrateMonomialsFaceReduction_(
    int c, int f, double factor, int k, Polynomial& integrals) const
{};


/* ******************************************************************
* Create surface integration object
****************************************************************** */
Polynomial ConvertPolynomialsToSurfacePolynomial(
    const AmanziGeometry::Point& xf, 
    const std::shared_ptr<SurfaceCoordinateSystem>& coordsys,
    const std::vector<const PolynomialBase*>& polys)
{
  int d = xf.dim();
  Polynomial product(d - 1, 0);
  product(0) = 1.0;

  for (int i = 0; i < polys.size(); ++ i) {
    Polynomial tmp(d, polys[i]->order(), polys[i]->ExpandCoefficients());
    tmp.set_origin(polys[i]->get_origin());
    tmp.ChangeCoordinates(xf, *coordsys->tau());  
    product *= tmp;
  }

  return product;
}

}  // namespace WhetStone
}  // namespace Amanzi

