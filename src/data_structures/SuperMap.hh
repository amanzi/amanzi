/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! SuperMap class provides a convenient way of creating and using SuperMapLumped
/*
  Amanzi uses SuperMapLumped in a few different ways that make its natural interface
  not that convenient.  SuperMapLumped also has a few limitations that simplify its
  design and implementaiton a great deal.

  This class is a Helper class in that it wraps a SuperMapLumped, providing a
  related interface that is better designed for users.

  It also enforces and mitigates the design limitations of SuperMapLumped itself.
*/

#ifndef AMANZI_OPERATORS_SUPER_MAP_WRAPPER_HH_
#define AMANZI_OPERATORS_SUPER_MAP_WRAPPER_HH_

#include <map>
#include "Teuchos_RCP.hpp"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "SuperMapLumped.hh"

class Epetra_Vector;

namespace Amanzi {

class CompositeVectorSpace;
class CompositeVector;
class TreeVectorSpace;
class TreeVector;

//namespace AmanziMesh {
//class Mesh;
//}


namespace Operators {

// wrapper class
class SuperMap {
 public:
  SuperMap(const std::vector<CompositeVectorSpace>& cvss);

  // map accessors
  // -- global map accessors
  Teuchos::RCP<const Epetra_Map> Map() const { return smap_->Map(); }
  Teuchos::RCP<const Epetra_Map> GhostedMap() const { return smap_->GhostedMap(); }

  // -- component map accessors
  Teuchos::RCP<const Epetra_BlockMap> ComponentMap(int block_num, const std::string& compname) const
  {
    return smap_->ComponentMap(block_info_.at(std::make_tuple(block_num, compname, 0)).first);
  }

  Teuchos::RCP<const Epetra_BlockMap>
  ComponentGhostedMap(int block_num, const std::string& compname) const
  {
    return smap_->ComponentGhostedMap(
      block_info_.at(std::make_tuple(block_num, compname, 0)).first);
  }

  // check if accessor is valid
  bool HasComponent(int block_num, const std::string& compname, int dof_num = 0) const
  {
    return block_info_.count(std::make_tuple(block_num, compname, dof_num)) != 0;
  }

  // index accessors
  const std::vector<int>& Indices(int block_num, const std::string& compname, int dof_num) const
  {
    auto bi = block_info_.find(std::make_tuple(block_num, compname, dof_num));
    if (bi == block_info_.end()) {
      Errors::Message msg;
      msg << "SuperMap does not have block component <" << block_num << "," << compname << ","
          << dof_num << ">";
      Exceptions::amanzi_throw(msg);
    }
    return smap_->Indices(bi->second.first, bi->second.second);
  }

  const std::vector<int>&
  GhostIndices(int block_num, const std::string& compname, int dof_num) const
  {
    auto bi = block_info_.find(std::make_tuple(block_num, compname, dof_num));
    if (bi == block_info_.end()) {
      Errors::Message msg;
      msg << "SuperMap does not have block component <" << block_num << "," << compname << ","
          << dof_num << ">";
      Exceptions::amanzi_throw(msg);
    }
    return smap_->GhostIndices(bi->second.first, bi->second.second);
  }

  // block indices.  This is an array of integers, length Map().MyLength(),
  // where each dof and component have a unique integer value.  The returned
  // int is the number of unique values, equal to
  // sum(NumDofs(comp) for comp in components), in this array.
  std::pair<int, Teuchos::RCP<std::vector<int>>> BlockIndices() const
  {
    return smap_->BlockIndices();
  }

 protected:
  std::unique_ptr<SuperMapLumped> smap_;
  std::map<std::tuple<int, std::string, int>, std::pair<std::string, int>> block_info_;
};


// Nonmember contructors/factories
Teuchos::RCP<SuperMap>
createSuperMap(const CompositeVectorSpace& cv);
Teuchos::RCP<SuperMap>
createSuperMap(const TreeVectorSpace& cv);

// Copy in/out
int
copyToSuperVector(const SuperMap& map,
                  const CompositeVector& bv,
                  Epetra_Vector& sv,
                  int block_num = 0);
int
copyFromSuperVector(const SuperMap& map,
                    const Epetra_Vector& sv,
                    CompositeVector& bv,
                    int block_num = 0);
int
addFromSuperVector(const SuperMap& map,
                   const Epetra_Vector& sv,
                   CompositeVector& bv,
                   int block_num = 0);

// Nonmember TreeVector to/from Super-vector
// -- simple schema version
int
copyToSuperVector(const SuperMap& map, const TreeVector& tv, Epetra_Vector& sv);
int
copyFromSuperVector(const SuperMap& map, const Epetra_Vector& sv, TreeVector& tv);
int
addFromSuperVector(const SuperMap& map, const Epetra_Vector& sv, TreeVector& tv);

} // namespace Operators
} // namespace Amanzi

#endif
