/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Multiphase PK

  Boundary function.
*/

#ifndef AMANZI_MULTIPHASE_BOUNDARY_FUNCTION_HH_
#define AMANZI_MULTIPHASE_BOUNDARY_FUNCTION_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "PK_DomainFunction.hh"

#include "MultiphaseDefs.hh"

namespace Amanzi {
namespace Multiphase {

class MultiphaseBoundaryFunction : public PK_DomainFunction {
 public:
  MultiphaseBoundaryFunction()
    : rainfall_(false),
      bc_name_("underfined"),
      component_id_(-1),
      component_name_("water"),
      component_phase_(-1) {};

  MultiphaseBoundaryFunction(const Teuchos::ParameterList& plist);
  void ComputeSubmodel(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // modifiers and access
  void set_bc_name(const std::string& name) { bc_name_ = name; }
  std::string get_bc_name() { return bc_name_; }

  void SetComponentId(const std::vector<std::string>& names);
  std::string component_name() { return component_name_; }
  int component_id() { return component_id_; }
  int component_phase() { return component_phase_; }

 private:
  bool rainfall_;
  std::string bc_name_;

  int component_id_;
  std::string component_name_;
  int component_phase_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
