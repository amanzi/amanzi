/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Mesh Functions

  This is a generalization of a mesh function: the function domain
  does not match a mesh and volume fraction are introduced.

  Known issues: if function domain is defined by a set of overlapping
  regions, the volume fractions are summed up which is correct only
  for special configurations.
*/

#ifndef AMANZI_MATERIAL_MESH_FUNCTION_HH_
#define AMANZI_MATERIAL_MESH_FUNCTION_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "MeshFunction.hh"

namespace Amanzi {
namespace Functions {

class MaterialMeshFunction : public MeshFunction {
 public:
  MaterialMeshFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : MeshFunction(mesh){};

  // Overload the AddSpec method to check uniqueness and sum-up
  // volume fractions.
  virtual void AddSpec(const Teuchos::RCP<Spec>& spec);

 protected:
  typedef std::map<AmanziMesh::Entity_ID, double> MaterialMesh;

  typedef std::pair<Teuchos::RCP<Spec>, Teuchos::RCP<MaterialMesh>> MaterialSpec;
  typedef std::vector<Teuchos::RCP<MaterialSpec>> MaterialSpecList;

  std::map<AmanziMesh::Entity_kind, Teuchos::RCP<MaterialSpecList>> material_specs_;
};

} // namespace Functions
} // namespace Amanzi

#endif
