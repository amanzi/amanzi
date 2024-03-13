/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

/*
 Shallow water PK

*/

#ifndef AMANZI_SHALLOW_WATER_BOUNDARY_FUNCTION_HH_
#define AMANZI_SHALLOW_WATER_BOUNDARY_FUNCTION_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "WhetStoneDefs.hh"

#include "PK_DomainFunction.hh"

namespace Amanzi {
namespace ShallowWater {

class ShallowWaterBoundaryFunction : public PK_DomainFunction {
 public:
  ShallowWaterBoundaryFunction() : bc_name_("undefined"){};
  ShallowWaterBoundaryFunction(const Teuchos::ParameterList& plist);

  // modifiers and access
  void set_bc_name(const std::string& name) { bc_name_ = name; }
  std::string get_bc_name() { return bc_name_; }

  void set_type(WhetStone::DOF_Type type) { type_ = type; }
  WhetStone::DOF_Type type() { return type_; }

  std::vector<double> bc_value(int f) { return value_[f]; }
  bool bc_find(int f) { return value_.find(f) != value_.end(); }

 private:
  std::string bc_name_;
  WhetStone::DOF_Type type_; // type of dofs related to this bc

  std::vector<std::string> regions_;
};

} // namespace ShallowWater
} // namespace Amanzi

#endif
