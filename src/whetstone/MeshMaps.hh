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
#include "MatrixObjects.hh"
#include "Tensor.hh"
#include "SpaceTimePolynomial.hh"
#include "VectorObjects.hh"

namespace Amanzi {
namespace WhetStone {

class Polynomial;

class MeshMaps {
 public:
  MeshMaps(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : mesh0_(mesh), mesh1_(mesh), d_(mesh1_->getSpaceDimension()){};

  MeshMaps(Teuchos::RCP<const AmanziMesh::Mesh> mesh0, Teuchos::RCP<const AmanziMesh::Mesh> mesh1)
    : mesh0_(mesh0), mesh1_(mesh1), d_(mesh1_->getSpaceDimension()){};

  virtual ~MeshMaps(){};

  // Maps
  // -- pseudo-velocity
  virtual void VelocityEdge(int e, VectorPolynomial& ve) const;
  virtual void VelocityFace(int f, VectorPolynomial& vf) const;
  virtual void VelocityCell(int c,
                            const std::vector<VectorPolynomial>& ve,
                            const std::vector<VectorPolynomial>& vf,
                            VectorPolynomial& vc) const = 0;

  // -- Nanson formula. Face deformation is defined completely by the
  //    deformation map in this formula: X = x + map(x)
  void
  NansonFormula(int f, const VectorSpaceTimePolynomial& map, VectorSpaceTimePolynomial& cn) const;

  // -- Jacobian
  void Jacobian(const VectorPolynomial& vc, MatrixPolynomial& J) const;

  // -- matrix of cofactors
  template <typename Matrix>
  void Cofactors(const Matrix& J, Matrix& C) const;

  // -- determinant
  template <typename Matrix, typename Vector>
  void Determinant(const Matrix& J, Vector& det) const;

  // Miscalleneous
  // -- projection ffrom reference coordinates (mesh0) to mesh1
  void ProjectPolynomial(int c, Polynomial& poly) const;

  // -- polynomial approximation of map x2 = F(x1)
  int LeastSquareFit(int order,
                     const AmanziMesh::Point_List& x1,
                     const AmanziMesh::Point_List& x2,
                     VectorPolynomial& u) const;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh0_; // initial mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh1_; // target mesh
  int d_;
};


/* ******************************************************************
* Calculation of matrix of cofactors.
****************************************************************** */
template <typename Matrix>
void
MeshMaps::Cofactors(const Matrix& J, Matrix& C) const
{
  // allocate memory for matrix of cofactors
  C.Reshape(d_, d_, d_, 0, false);

  // calculate cofactors
  if (d_ == 2) {
    C(1, 1) = J(0, 0);
    C(1, 0) = J(0, 1);
    C(1, 0) *= -1.0;

    C(0, 0) = J(1, 1);
    C(0, 1) = J(1, 0);
    C(0, 1) *= -1.0;
  } else if (d_ == 3) {
    C(0, 0) = J(1, 1) * J(2, 2) - J(2, 1) * J(1, 2);
    C(1, 0) = J(2, 1) * J(0, 2) - J(0, 1) * J(2, 2);
    C(2, 0) = J(0, 1) * J(1, 2) - J(1, 1) * J(0, 2);

    C(0, 1) = J(2, 0) * J(1, 2) - J(1, 0) * J(2, 2);
    C(1, 1) = J(0, 0) * J(2, 2) - J(2, 0) * J(0, 2);
    C(2, 1) = J(1, 0) * J(0, 2) - J(0, 0) * J(1, 2);

    C(0, 2) = J(1, 0) * J(2, 1) - J(2, 0) * J(1, 1);
    C(1, 2) = J(2, 0) * J(0, 1) - J(0, 0) * J(2, 1);
    C(2, 2) = J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1);
  }
}


/* ******************************************************************
* Calculate detminant of a matrix.
****************************************************************** */
template <typename Matrix, typename Poly>
void
MeshMaps::Determinant(const Matrix& J, Poly& det) const
{
  if (d_ == 2) {
    det = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
  } else if (d_ == 3) {
    det = J(0, 0) * J(1, 1) * J(2, 2) + J(2, 0) * J(0, 1) * J(1, 2) + J(1, 0) * J(2, 1) * J(0, 2) -
          J(2, 0) * J(1, 1) * J(0, 2) - J(1, 0) * J(0, 1) * J(2, 2) - J(0, 0) * J(2, 1) * J(1, 2);
  }
}

} // namespace WhetStone
} // namespace Amanzi

#endif
