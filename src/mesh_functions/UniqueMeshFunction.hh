/*
  Mesh Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Function applied to a mesh component with at most one function
  application per entity.
*/

#ifndef AMANZI_UNIQUE_MESH_FUNCTION_HH_
#define AMANZI_UNIQUE_MESH_FUNCTION_HH_

#include <vector>
#include <set>
#include <string>

#include "Teuchos_RCP.hpp"

#include "MeshFunction.hh"

namespace Amanzi {
namespace Functions {

class UniqueMeshFunction : public MeshFunction {
 public:
  // Constructor
  UniqueMeshFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : MeshFunction(mesh){};

  // Overload the AddSpec method to check uniqueness.
  virtual void AddSpec(const Teuchos::RCP<Spec>& spec);

 protected:
  typedef std::set<AmanziMesh::Entity_ID> MeshIDs;
  typedef std::pair<Teuchos::RCP<Spec>, Teuchos::RCP<MeshIDs>> UniqueSpec;
  typedef std::vector<Teuchos::RCP<UniqueSpec>> UniqueSpecList;

  std::map<AmanziMesh::Entity_kind, Teuchos::RCP<UniqueSpecList>> unique_specs_;
};

} //namespace Functions
} //namespace Amanzi

#endif
