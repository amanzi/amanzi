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

  A polynomial binded with a mesh object. This struct allows us
  to verify identity of a polynomial used by multiple classes.
*/

#ifndef AMANZI_WHETSTONE_POLYNOMIAL_ON_MESH_HH_
#define AMANZI_WHETSTONE_POLYNOMIAL_ON_MESH_HH_

#include "Teuchos_RCP.hpp"

#include "Point.hh"

#include "Polynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class PolynomialOnMesh {
 public:
  PolynomialOnMesh()
    : kind_((Entity_kind)WhetStone::CELL), id_(-1) {};

  Polynomial& poly() { return poly_; }
  const Polynomial& poly() const { return poly_; }

  void set_kind(Entity_kind kind) { kind_ = kind; }
  const Entity_kind& get_kind() const { return kind_; }

  const Entity_ID& get_id() const { return id_; }
  void set_id(Entity_ID id) { id_ = id; }

 private:
  Polynomial poly_;
  Entity_kind kind_; // topological binding of polynomial
  Entity_ID id_;     // numerical id of topological entity
};

} // namespace WhetStone
} // namespace Amanzi

#endif
