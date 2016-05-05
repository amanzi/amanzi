/*
  Mesh Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov

  Function applied to a mesh component defined on a material region
  with at most one function application per entity.
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
  MaterialMeshFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
      : MeshFunction(mesh) {};

  // Overload the AddSpec method to check uniqueness and sum-up 
  // volume fractions.
  virtual void AddSpec(const Teuchos::RCP<Spec>& spec);

 protected:
  typedef std::map<AmanziMesh::Entity_ID, double> MaterialMesh;

  typedef std::pair<Teuchos::RCP<Spec>, Teuchos::RCP<MaterialMesh> > UniqueSpec;
  typedef std::vector<Teuchos::RCP<UniqueSpec> > UniqueSpecList;

  std::map<AmanziMesh::Entity_kind, Teuchos::RCP<UniqueSpecList> > unique_specs_;
};

}  // namespace Functions
}  // namespace Amanzi

#endif
