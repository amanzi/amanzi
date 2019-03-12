/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

//! SuperMapWrapper class provides a convenient way of creating and using SuperMap

/*
  Amanzi uses SuperMap in a few different ways that make its natural interface
  not that convenient.  SuperMap also has a few limitations that simplify its
  design and implementaiton a great deal.

  This class is a Helper class in that it wraps a SuperMap, providing a
  related interface that is better designed for users.

  It also enforces and mitigates the design limitations of SuperMap itself.
*/

#ifndef AMANZI_OPERATORS_SUPER_MAP_WRAPPER_HH_
#define AMANZI_OPERATORS_SUPER_MAP_WRAPPER_HH_


namespace Amanzi {

class CompositeVectorSpace;
class TreeVectorSpace;

namespace Operators {


// Nonmember contructors/factories
Teuchos::RCP<SuperMapWrapper> createSuperMap(const CompositeVectorSpace& cv);
Teuchos::RCP<SuperMapWrapper> createSuperMap(const TreeVectorSpace& cv);

// Nonmember helper function
std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> >
getMaps(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_kind location);


// wrapper class
class SuperMapWrapper {
 public:
  SuperMapWrapper(const std::vector<CompositeVectorSpace>& cvss);

  // map accessors
  // -- global map accessors
  Teuchos::RCP<const Epetra_Map> Map() const { return smap_->Map(); }
  Teuchos::RCP<const Epetra_Map> GhostedMap() const { return smap_->GhostedMap(); }

  // -- component map accessors
  Teuchos::RCP<const Epetra_BlockMap>
  ComponentMap(int block_num, const std::string& compname, int dof_num);

  Teuchos::RCP<const Epetra_BlockMap>
  ComponentGhostedMap(int block_num, const std::string& compname, int dof_num);
  
  // index accessors
  const std::vector<int>&
  Indices(int block_num, const std::string& compname, int dof_num) const;
  const std::vector<int>&
  GhostIndices(int block_num, const std::string& compname, int dofnum) const;

 protected:
  std::unique_ptr<SuperMap> smap_;


};


} // namespace
} // namespace



#endif
