/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

//! SuperMap class for blocking same-type maps into a single, common map.

/*!
  Takes non-contiguous data structure spaces (CompositeVector, TreeVector)
  and converts them into a single map.

  DESIGN FLAW: Currently this assumes that component names are unique, and if
  two components share the same name, they share the same map.  This is
  obviously wrong when multple meshes are involved -- for instance a TV of
  surface + subsurface, both with "cell" components, would break miserably.

  That said, this is a simple and well-posed class.  It stays, but users
  should prefer to use the wrapper class, SuperMapWrapper, which deals with a
  few shortcomings of this class.
*/

#ifndef AMANZI_OPERATORS_SUPER_MAP_HH_
#define AMANZI_OPERATORS_SUPER_MAP_HH_

#include "Teuchos_RCP.hpp"
#include "AmanziComm.hh"

#include "dbc.hh"
#include "Mesh.hh"


namespace Amanzi {
namespace Operators {

class SuperMap {
 public:
  // Constructor
  SuperMap(const Comm_ptr_type& comm,
           const std::vector<std::string>& compnames,
           const std::vector<int>& dofnums,
           const std::vector<Teuchos::RCP<const Epetra_BlockMap> >& maps,
           const std::vector<Teuchos::RCP<const Epetra_BlockMap> >& ghost_maps);

  SuperMap(const SuperMap& other);  
  virtual ~SuperMap() = default;
  
  // map accessors
  Teuchos::RCP<const Epetra_Map> Map() const { return map_; }
  Teuchos::RCP<const Epetra_Map> GhostedMap() const { return ghosted_map_; }

  // index accessors
  const std::vector<int>& Indices(const std::string& compname, int dofnum) const;
  const std::vector<int>& GhostIndices(const std::string& compname, int dofnum) const;

  // block indices.  This is an array of integers, length Map().MyLength(),
  // where each dof and component have a unique integer value.  The returned
  // int is the number of unique values, equal to
  // sum(NumDofs(comp) for comp in components), in this array.
  std::pair<int, Teuchos::RCP<std::vector<int> > > BlockIndices() const;

#ifdef SUPERMAP_TESTING
 public:
#else
 protected:
#endif
  
  // meta-data accessors
  bool HasComponent(const std::string& compname) const;
  int Offset(const std::string& compname) const { return offsets_.at(compname); }
  int GhostedOffset(const std::string& compname) const { return ghosted_offsets_.at(compname); }
  int NumOwnedElements(const std::string& compname) const { return counts_.at(compname); }
  int NumUsedElements(const std::string& compname) const {
    return counts_.at(compname) + ghosted_counts_.at(compname); }
  int NumDofs(const std::string& compname) const { return num_dofs_.at(compname); }

  // iterate over compnames
  typedef std::vector<std::string>::const_iterator name_iterator;
  name_iterator begin() const { return compnames_.begin(); }
  name_iterator end() const { return compnames_.end(); }
  unsigned int size() const { return compnames_.size(); }

 protected:
  virtual const std::vector<int>& CreateIndices_(const std::string& compname, int dofnum, bool ghosted) const;

 protected:
  std::vector<std::string> compnames_;
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


} // namespace Operators
} // namespace Amanzi

#endif
