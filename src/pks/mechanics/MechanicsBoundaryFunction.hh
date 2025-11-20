/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
 Shallow water PK

*/

#ifndef AMANZI_MECHANICS_BOUNDARY_FUNCTION_HH_
#define AMANZI_MECHANICS_BOUNDARY_FUNCTION_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "WhetStoneDefs.hh"

#include "PK_DomainFunction.hh"

namespace Amanzi {
namespace Mechanics {

class MechanicsBoundaryFunction : public PK_DomainFunction {
 public:
  MechanicsBoundaryFunction()
    : bc_name_("undefined") {};
  MechanicsBoundaryFunction(const Teuchos::ParameterList& plist);

  // modifiers and access
  void set_bc_name(const std::string& name) { bc_name_ = name; }
  std::string get_bc_name() { return bc_name_; }

  void set_type(WhetStone::DOF_Type type) { type_ = type; }
  WhetStone::DOF_Type type() { return type_; }

  void set_kind(AmanziMesh::Entity_kind kind) { kind_ = kind; }
  AmanziMesh::Entity_kind kind() { return kind_; }

  std::string plane_strain_direction() { return plane_strain_direction_; }

 private:
  std::string bc_name_;
  WhetStone::DOF_Type type_; // type of dofs related to this bc
  AmanziMesh::Entity_kind kind_;

  std::string plane_strain_direction_;
};

} // namespace Mechanics
} // namespace Amanzi

#endif
