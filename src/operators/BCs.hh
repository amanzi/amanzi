/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#ifndef AMANZI_OPERATORS_BC_HH_
#define AMANZI_OPERATORS_BC_HH_

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "WhetStoneDefs.hh"

#include "OperatorDefs.hh"

namespace Amanzi {

namespace AmanziGeometry {
class Point;
}

namespace Operators {

/* *******************************************************************
* Elliptic equation E(u) = f. Three types of boundary conditions are
* supported by this class:
*   [Dirichlet]                  u = u0
*   [Neumann]     -K(u) grad u . n = g0
*   [Mixed] -K(u) grad u . n - c u = g1
*
* The right-hand side data (u0, g0, g1) must be placed in array
* bc_value that has a proper size (see below). The type of BC
* must be indicated in integer array bc_model using constants
* defined in file OperatorDefs.hh. Arrays bc_value and bc_model
* must have the same size and contain ghost degrees of freedom.
*
* The coefficent c must be placed in array bc_mixed. This array
* can be empty; otherwise, its size must match that of bc_value.
*
* All three arrays are associated with degrees of freedom selected
* for a problem discretization, see class Operators for more detail.
* For example, for the nodal discretization of elliptic equation,
* the dimension of the arrays equals to the total number of nodes
* on a processor, including the ghost nodes.
*
* NOTE. Arrays bc_value and bc_model may be empty when homogeneous
*   Neumann boundary conditions are imposed on the domain boundary.
*
* NOTE. Suffient conditions for solution non-negativity are
*   g0 <= 0, g1 <= 0 and c >=0.
*
* NOTE. All data in input arrays are given with respect to exterior
*   normal vector. Implementation of boundary conditions should take
*   into account that actual mesh normal may be oriented arbitrarily.
*
* **********************
*
* Diffusion-advection equation E(u) + A(u) = f. Four types of boundary
* conditions are supported:
*   [Dirichlet]                         u = u0
*   [Neumann]            -K(u) grad u . n = g0
*   [Mixed]        -K(u) grad u . n - c u = g1
*   [Total flux] -(K(u) grad u - v c) . n = g2
*
* Here v is the advective velocity. For the diffusion-advection
* operator, we may impose boundary conditions that make sence for
* diffusion but not appropriate for advection. To void creation of two
* sets of boundary conditions, the total flux condition can be used.
* Only the leading operator, typically diffusion, can set up this BC.
* The other operators will remove all boundary contributions to the
* matrix and right-hand side when the total flux condition is specified.
*
* **********************
*
* Advection equation A(u) = f. One type of boundary condition is
* supported by this class:
*   [Dirichlet]          u = u0
*
* The data u0 can be included in a weak formulation in two different
* ways. In the integration by parts formulations, boundary integrals
* are elliminated from the left-hand side and the right-hand side is
* modified, so this boundary conditions behaves similar to the
* diffusion problem. For other weak formulations, a weak form of this
* boundary condition is added to the system:
*
*   (A(u) - f, w) + (u - u0, w) = 0.
*
* In the second approach, array bc_model should use TYPE2 boundary
* condition, see OperatorDefs.hh for the full name.
*
* NOTE. When no contibition from boundary should occur, array bc_model
* must contain integer value OPERATOR_BC_REMOVE.
******************************************************************* */

class BCs {
 public:
  // KIND defines location of DOFs on a mesh:
  // -- available mesh entities are node, edge, face, and cell
  //
  // TYPE provides additional information, see OperatorDefs.hh for available options.
  // In short, is specifies geometric, algebraic or any other information:
  // -- scalar is the simplest DOF, it is just a number (example: mean pressure)
  // -- point is a vector DOF which has mesh dimension (example: fluid velocity)
  // -- vector is a general vector DOF (example: moments of pressure)
  // -- normal-component is a geometric DOF (example: normal component of fluid velocity)
  BCs(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
      AmanziMesh::Entity_kind kind,
      WhetStone::DOF_Type type)
    : kind_(kind), type_(type), mesh_(mesh){};
  ~BCs(){};

  // non-const access
  Teuchos::RCP<const AmanziMesh::Mesh> mesh() const { return mesh_; }
  AmanziMesh::Entity_kind kind() const { return kind_; }
  WhetStone::DOF_Type type() const { return type_; }

  std::vector<int>& bc_model()
  {
    if (bc_model_.size() == 0) {
      int nent = mesh_->getNumEntities(kind_, AmanziMesh::Parallel_type::ALL);
      bc_model_.resize(nent, Operators::OPERATOR_BC_NONE);
    }
    return bc_model_;
  }

  std::vector<double>& bc_value()
  {
    if (bc_value_.size() == 0) {
      int nent = mesh_->getNumEntities(kind_, AmanziMesh::Parallel_type::ALL);
      bc_value_.resize(nent, 0.0);
    }
    return bc_value_;
  }

  std::vector<double>& bc_mixed()
  {
    if (bc_mixed_.size() == 0) {
      int nent = mesh_->getNumEntities(kind_, AmanziMesh::Parallel_type::ALL);
      bc_mixed_.resize(nent, 0.0);
    }
    return bc_mixed_;
  }

  std::vector<AmanziGeometry::Point>& bc_value_point()
  {
    if (bc_value_point_.size() == 0) {
      AmanziGeometry::Point p(mesh_->getSpaceDimension());
      int nent = mesh_->getNumEntities(kind_, AmanziMesh::Parallel_type::ALL);
      bc_value_point_.resize(nent, p);
    }
    return bc_value_point_;
  }

  std::vector<std::vector<double>>& bc_value_vector(int n = 1)
  {
    if (bc_value_vector_.size() == 0) {
      int nent = mesh_->getNumEntities(kind_, AmanziMesh::Parallel_type::ALL);
      bc_value_vector_.resize(nent);

      for (int i = 0; i < nent; ++i) { bc_value_vector_[i].resize(n); }
    }
    return bc_value_vector_;
  }

  // const access
  const std::vector<int>& bc_model() const { return bc_model_; }
  const std::vector<double>& bc_value() const { return bc_value_; }
  const std::vector<double>& bc_mixed() const { return bc_mixed_; }
  const std::vector<AmanziGeometry::Point>& bc_value_point() const { return bc_value_point_; }
  const std::vector<std::vector<double>>& bc_value_vector() const { return bc_value_vector_; }

 private:
  AmanziMesh::Entity_kind kind_;
  WhetStone::DOF_Type type_;

  std::vector<int> bc_model_;
  std::vector<double> bc_value_;
  std::vector<double> bc_mixed_;

  std::vector<std::vector<double>> bc_value_vector_;
  std::vector<AmanziGeometry::Point> bc_value_point_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

} // namespace Operators
} // namespace Amanzi

#endif
