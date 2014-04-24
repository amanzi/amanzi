/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Function applied to a mesh component with at most one function application per
entity.

------------------------------------------------------------------------- */

#ifndef AMANZI_UNIQUE_MESH_FUNCTION_HH_
#define AMANZI_UNIQUE_MESH_FUNCTION_HH_

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "mesh_function.hh"

namespace Amanzi {
namespace Functions {

class UniqueMeshFunction : public MeshFunction {

public:
  // Constructor
  UniqueMeshFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    MeshFunction(mesh) {}

  // Overload the AddSpec method to check uniqueness.
  virtual void AddSpec(const Teuchos::RCP<Spec>& spec);

protected:
  typedef std::set<AmanziMesh::Entity_ID> SpecIDs;
  typedef std::pair<Teuchos::RCP<Spec>, Teuchos::RCP<SpecIDs> > SpecAndIDs;
  typedef std::vector<Teuchos::RCP<SpecAndIDs> > SpecAndIDsList;

  std::map<AmanziMesh::Entity_kind, Teuchos::RCP<SpecAndIDsList> > specs_and_ids_;

};

} //namespace
} //namespace

#endif
