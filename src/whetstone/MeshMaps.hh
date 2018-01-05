/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for maps between mesh objects located on different 
  meshes, e.g. two states (mesh0 and mesh1) of a deformable mesh. 
  For calculating time dependent maps, we connect the two states 
  using the linearized map: x + t (F(x) - x) where F(x) is a 
  steady-state map from mesh0 to mesh1 and time t is between 0 
  and 1. Thus, the linearized velocity is v = F(x) - x.
*/

#ifndef AMANZI_WHETSTONE_MESH_MAPS_HH_
#define AMANZI_WHETSTONE_MESH_MAPS_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "Tensor.hh"
#include "VectorPolynomial.hh"
#include "WhetStone_typedefs.hh"

namespace Amanzi {
namespace WhetStone {

class Polynomial;

class MeshMaps { 
 public:
  MeshMaps(Teuchos::RCP<const AmanziMesh::Mesh> mesh) 
    : mesh1_(mesh),
      mesh0_(mesh),
      d_(mesh1_->space_dimension()) {};

  MeshMaps(Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
           Teuchos::RCP<const AmanziMesh::Mesh> mesh1) 
    : mesh1_(mesh1),
      mesh0_(mesh0),
      d_(mesh1_->space_dimension()) {};

  ~MeshMaps() {};

  // Maps
  // -- pseudo-velocity
  virtual void VelocityFace(int f, VectorPolynomial& v) const;
  virtual void VelocityCell(int c, const std::vector<VectorPolynomial>& vf,
                            VectorPolynomial& vc) const = 0;

  // -- Nanson formula
  virtual void NansonFormula(int f, double t, const VectorPolynomial& vf,
                             VectorPolynomial& cn) const = 0;

  // -- Jacobian
  void Jacobian(const VectorPolynomial& vc, MatrixPolynomial& J) const;

  // -- matrix of cofactors
  void Cofactors(double t, const MatrixPolynomial& J, MatrixPolynomial& C) const;

  // -- determinant
  void Determinant(double t, const MatrixPolynomial& J, Polynomial& det) const;

  // Miscalleneous
  // -- polynomial approximation of map x2 = F(x1)
  int LeastSquareFit(int order,
                     const std::vector<AmanziGeometry::Point>& x1, 
                     const std::vector<AmanziGeometry::Point>& x2,
                     VectorPolynomial& u) const;

  // extension of mesh interface
  AmanziGeometry::Point cell_geometric_center(int id, int c) const;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh0_;  // initial mesh 
  Teuchos::RCP<const AmanziMesh::Mesh> mesh1_;  // target mesh
  int d_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

