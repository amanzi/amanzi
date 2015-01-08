/*
Author: Ethan Coon (ecoon@lanl.gov)

Takes non-contiguous data structure spaces (CompositeVector, TreeVector) and
converts them into a single map.

DESIGN FLAW: Currently this assumes that component names are unique, and if
two components share the same name, they share thee same map.  This is
obviously wrong when multple meshes are involved -- for instance a TV of
surface + subsurface, both with "cell" components, would break miserably.

*/

#ifndef AMANZI_OPERATORS_SUPER_MAP_HH_
#define AMANZI_OPERATORS_SUPER_MAP_HH_

#include "Teuchos_RCP.hpp"

namespace Amanzi {

class CompositeVectorSpace;
class TreeVectorSpace;

namespace Operators {

class SuperMap {

 public:

  // Constructor
  SuperMap(const Epetra_MpiComm& comm,
           const std::vector<std::string>& compnames,
           const std::vector<int>& dofnums,
           const std::vector<Teuchos::RCP<const Epetra_Map> >& maps,
           const std::vector<Teuchos::RCP<const Epetra_Map> >& ghost_maps);

  // meta-data accessors
  int Offset(std::string compname) const { return offsets_.at(compname); }
  int GhostedOffset(std::string compname) const { return ghosted_offsets_.at(compname); }
  int NumOwnedElements(std::string compname) const { return counts_.at(compname); }
  int NumUsedElements(std::string compname) const { return counts_.at(compname) + ghosted_counts_.at(compname); }
  int NumDofs(std::string compname) const { return num_dofs_.at(compname); }

  // map accessors
  Teuchos::RCP<const Epetra_Map> Map() const { return map_; }
  Teuchos::RCP<const Epetra_Map> GhostedMap() const { return ghosted_map_; }

  // index accessors
  const std::vector<int> Indices(std::string compname, int dofnum) const;
  const std::vector<int> GhostIndices(std::string compname, int dofnum) const;

 protected:
    const std::vector<int> CreateIndices_(std::string compname, int dofnum, bool ghosted) const;

 protected:
  std::map<std::string,int> offsets_;
  std::map<std::string,int> num_dofs_;
  std::map<std::string,int> counts_;
  std::map<std::string,int> ghosted_offsets_;
  std::map<std::string,int> ghosted_counts_;

  mutable std::map<std::string, std::map<int, std::vector<int> > > indices_;
  mutable std::map<std::string, std::map<int, std::vector<int> > > ghosted_indices_;

  Teuchos::RCP<Epetra_Map> map_;
  Teuchos::RCP<Epetra_Map> ghosted_map_;
};


// Nonmember helper function
std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> >
getMaps_(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_kind location);

// Nonmember contructors/factories
Teuchos::RCP<SuperMap>
createSuperMap(const CompositeVectorSpace& cv);

Teuchos::RCP<SuperMap>
createSuperMap(const TreeVectorSpace& cv);



} // namespace Operators
} // namespace Amanzi

#endif
