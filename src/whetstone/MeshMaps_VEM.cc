/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Maps between mesh objects located on different meshes, e.g. two states
  of a deformable mesh: virtual element implementation.
*/

#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D_CrouzeixRaviart.hh"
#include "MeshMaps_VEM.hh"
#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Calculate mesh velocity in cell c: new algorithm.
* NOTE: second mesh is not used, so does it belong here?
****************************************************************** */
void
MeshMaps_VEM::VelocityCell(int c,
                           const std::vector<VectorPolynomial>& ve,
                           const std::vector<VectorPolynomial>& vf,
                           VectorPolynomial& vc) const
{
  Teuchos::ParameterList plist;
  plist.set<std::string>("method", method_).set<int>("method order", order_);
  auto mfd = BilinearFormFactory::Create(plist, mesh0_);

  vc.resize(d_);

  if (projector_ == "least square") {
    LeastSquareProjector_Cell_(order_, c, vf, vc);
  } else {
    for (int i = 0; i < d_; ++i) {
      std::vector<Polynomial> vvf, vve;
      for (int n = 0; n < vf.size(); ++n) { vvf.push_back(vf[n][i]); }

      for (int n = 0; n < ve.size(); ++n) { vve.push_back(ve[n][i]); }

      if (projector_ == "H1") {
        mfd->H1Cell(c, vve, vvf, NULL, vc[i]);
      } else if (projector_ == "L2") {
        mfd->L2Cell(c, vve, vvf, NULL, vc[i]);
      }
    }
  }
}


/* ******************************************************************
* Calculate mesh velocity on face f.
****************************************************************** */
void
MeshMaps_VEM::VelocityFace(int f, VectorPolynomial& vf) const
{
  if (d_ == 2) {
    MeshMaps::VelocityFace(f, vf);
  } else {
    
    auto [edges, dirs] = mesh0_->getFaceEdgesAndDirections(f);
    int nedges = edges.size();

    Teuchos::ParameterList plist;
    plist.set<std::string>("method", method_).set<int>("method order", order_);
    auto mfd = BilinearFormFactory::Create(plist, mesh0_);
    mfd->set_order(order_);

    vf.resize(d_);
    for (int i = 0; i < d_; ++i) {
      VectorPolynomial v;
      std::vector<Polynomial> ve;

      for (int n = 0; n < nedges; ++n) {
        int e = edges[n];
        MeshMaps::VelocityEdge(e, v);
        ve.push_back(v[i]);
      }

      if (projector_ == "L2")
        mfd->L2Face(f, ve, NULL, vf[i]);
      else if (projector_ == "H1")
        mfd->H1Face(f, ve, NULL, vf[i]);
    }
  }
}


/* ******************************************************************
* Calculate mesh velocity in cell c: old algorithm
****************************************************************** */
void
MeshMaps_VEM::LeastSquareProjector_Cell_(int order,
                                         int c,
                                         const std::vector<VectorPolynomial>& vf,
                                         VectorPolynomial& vc) const
{
  AMANZI_ASSERT(order == 1 || order == 2);

  vc.Reshape(d_, d_, order);

  AmanziGeometry::Point px1, px2;
  AmanziMesh::Point_List x1, x2;

  auto nodes = mesh0_->getCellNodes(c);
  int nnodes = nodes.size();

  for (int n = 0; n < nnodes; ++n) {
    px1 = mesh0_->getNodeCoordinate(nodes[n]);
    x1.push_back(px1);

    px2 = mesh1_->getNodeCoordinate(nodes[n]);
    x2.push_back(px2 - px1);
  }

  // FIXME
  if (order > 1) {
    const auto& faces = mesh0_->getCellFaces(c);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; ++n) {
      const auto& xf = mesh0_->getFaceCentroid(faces[n]);
      x1.push_back(xf);

      for (int i = 0; i < d_; ++i) { px2[i] = vf[n][i].Value(xf); }
      x2.push_back(px2);
    }
  }

  // calculate velocity u(X) = F(X) - X
  LeastSquareFit(order, x1, x2, vc);
}


/* ******************************************************************
* Calculate mesh velocity in cell c: old algorithm
****************************************************************** */
void
MeshMaps_VEM::ParseInputParameters_(const Teuchos::ParameterList& plist)
{
  method_ = plist.get<std::string>("method");
  order_ = plist.get<int>("method order");
  projector_ = plist.get<std::string>("projector");
}

} // namespace WhetStone
} // namespace Amanzi
