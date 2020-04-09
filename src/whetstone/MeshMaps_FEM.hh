/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Maps between mesh objects located on different meshes, e.g. two
  states of a deformable mesh: finite element implementation.
  Note that all points and polynomials use standard reference
  coordinates.
*/

#ifndef AMANZI_WHETSTONE_MESH_MAPS_FEM_HH_
#define AMANZI_WHETSTONE_MESH_MAPS_FEM_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "MeshMaps.hh"
#include "Polynomial.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class MeshMaps_FEM : public MeshMaps {
 public:
  MeshMaps_FEM(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : MeshMaps(mesh){};
  MeshMaps_FEM(Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
               Teuchos::RCP<const AmanziMesh::Mesh> mesh1)
    : MeshMaps(mesh0, mesh1){};
  ~MeshMaps_FEM(){};

  // remap pseudo velocity
  virtual void VelocityCell(int c, const std::vector<VectorPolynomial>& vf,
                            VectorPolynomial& vc) const override;

 private:
  void JacobianCellValue_(int c, double t, const AmanziGeometry::Point& x,
                          Tensor& J) const;

  Tensor JacobianValueInternal_(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                                int c, const AmanziGeometry::Point& xref) const;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
