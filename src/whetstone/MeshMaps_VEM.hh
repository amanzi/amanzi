/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Maps between mesh objects located on different meshes, e.g. two 
  states of a deformable mesh: virtual element implementation.
*/

#ifndef AMANZI_WHETSTONE_MESH_MAPS_VEM_HH_
#define AMANZI_WHETSTONE_MESH_MAPS_VEM_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "MeshMaps.hh"
#include "Polynomial.hh"
#include "Tensor.hh"
#include "WhetStone_typedefs.hh"

namespace Amanzi {
namespace WhetStone {

class MeshMaps_VEM : public MeshMaps { 
 public:
  MeshMaps_VEM(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : MeshMaps(mesh) {};
  MeshMaps_VEM(Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
               Teuchos::RCP<const AmanziMesh::Mesh> mesh1) : MeshMaps(mesh0, mesh1) {};
  ~MeshMaps_VEM() {};

  // Maps
  // -- pseudo-velocity in cell c
  virtual void VelocityCell(int c, VectorPolynomial& v) const override;
  // -- pseudo-velocity on face f
  virtual void VelocityFace(int f, VectorPolynomial& v) const override;

  // Jacobian
  // -- determinant of Jacobian
  virtual void JacobianDet(int c, double t, const std::vector<VectorPolynomial>& vf,
                           Polynomial& vc) const;

  // -- Jacobian value at point x
  virtual void JacobianCellValue(int c,
                                 double t, const AmanziGeometry::Point& x,
                                 Tensor& J) const;
  virtual void JacobianFaceValue(int f, const VectorPolynomial& v,
                                 const AmanziGeometry::Point& x,
                                 Tensor& J) const;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

