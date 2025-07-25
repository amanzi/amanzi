/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Operators

*/

#ifndef AMANZI_OPERATORS_BC_FACTORY_HH_
#define AMANZI_OPERATORS_BC_FACTORY_HH_

#include "BCs.hh"

namespace Amanzi {
namespace Operators {

class BCs_Factory {
 public:
  BCs_Factory() {};

  Teuchos::RCP<const AmanziMesh::Mesh> mesh() const { return mesh_; }
  void set_mesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) { mesh_ = mesh; }
  void set_kind(AmanziMesh::Entity_kind kind) { kind_ = kind; }
  void set_type(WhetStone::DOF_Type type) { type_ = type; }

  Teuchos::RCP<BCs> Create() const
  {
    auto bc = Teuchos::rcp(new BCs(mesh_, kind_, type_));
    // these are called to force instantiation
    bc->bc_model();
    bc->bc_value();
    return bc;
  }

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  AmanziMesh::Entity_kind kind_;
  WhetStone::DOF_Type type_;
};

} // namespace Operators
} // namespace Amanzi


#endif
