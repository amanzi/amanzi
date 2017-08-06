/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Maps between mesh objects located possibly on different meshes.
*/

#ifndef AMANZI_WHETSTONE_MESH_MAPS_HH_
#define AMANZI_WHETSTONE_MESH_MAPS_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "Polynomial.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"
#include "WhetStone_typedefs.hh"

namespace Amanzi {
namespace WhetStone {

class MeshMaps { 
 public:
  MeshMaps(Teuchos::RCP<const AmanziMesh::Mesh> mesh) 
    : mesh1_(mesh),
      mesh0_(mesh),
      d_(mesh1_->space_dimension()),
      method_(WHETSTONE_METHOD_VEM) {};

  MeshMaps(Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
           Teuchos::RCP<const AmanziMesh::Mesh> mesh1) 
    : mesh1_(mesh1),
      mesh0_(mesh0),
      d_(mesh1_->space_dimension()),
      method_(WHETSTONE_METHOD_VEM) {};

  ~MeshMaps() {};

  // support of advection schemes
  void FaceVelocity(int c, int f, std::vector<Polynomial>& v) const;
  Tensor FaceJacobian(int c, int f, const std::vector<Polynomial>& v,
                      const AmanziGeometry::Point& x) const;

  // finite elements use three meshes: reference, initial and target
  AmanziGeometry::Point FEM_Map(int c, const AmanziGeometry::Point& xref) const;
  Tensor FEM_Jacobian(int c, const AmanziGeometry::Point& xref) const;
  void FEM_Jacobian(int c, Polynomial& jac) const;

  // miscalleneous
  void set_method(int method) { method_ = method; }

  // polynomial approximation of map x2 = F(x1)
  int LeastSquareFit(int order,
                     const std::vector<AmanziGeometry::Point>& x1, 
                     const std::vector<AmanziGeometry::Point>& x2,
                     std::vector<AmanziGeometry::Point>& u) const;

 private:
  // finite element maps 
  Tensor FEM_JacobianInternal_(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                               int c, const AmanziGeometry::Point& xref) const;
  Tensor FEM_FaceJacobian_(int c, int f, const std::vector<Polynomial>& v,
                           const AmanziGeometry::Point& xref) const;

  // virtual element maps 
  void VEM_FaceVelocity_(int c, int f, std::vector<Polynomial>& v) const;
  Tensor VEM_FaceJacobian_(int c, int f, const std::vector<Polynomial>& v,
                           const AmanziGeometry::Point& x) const;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh1_;  // target mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh0_;  // initial mesh 
  int d_;
  int method_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

