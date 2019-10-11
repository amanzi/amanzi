/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  A polynomial binded with a mesh object. This struct allows us
  to verify identity of a polynomial used by multiple classes.
*/

#ifndef AMANZI_WHETSTONE_POLYNOMIAL_ON_MESH_HH_
#define AMANZI_WHETSTONE_POLYNOMIAL_ON_MESH_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "Polynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

struct PolynomialOnMesh {
 public:
  PolynomialOnMesh() : id_(-1), kind_((Entity_kind)WhetStone::CELL){};

  Polynomial& poly() { return poly_; }
  const Polynomial& poly() const { return poly_; }

  const Entity_kind& kind() const { return kind_; }
  const Entity_ID& id() const { return id_; }

  void set_kind(Entity_kind kind) { kind_ = kind; }
  void set_id(Entity_ID id) { id_ = id; }

 private:
  Polynomial poly_;
  Entity_kind kind_; // topological binding of polynomial
  Entity_ID id_;     // numerical id of topological entity
};

} // namespace WhetStone
} // namespace Amanzi

#endif
