/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

//! SuperMapLumped class for blocking same-type maps into a single, common map.

/*!
  Takes non-contiguous data structure spaces (CompositeVector, TreeVector)
  and converts them into a single map.

  DESIGN FLAW: Currently this assumes that component names are unique, and if
  two components share the same name, they share the same map.  This is
  obviously wrong when multple meshes are involved -- for instance a TV of
  surface + subsurface, both with "cell" components, would break miserably.

  That said, this is a simple and well-posed class.  It stays, but users
  should prefer to use the wrapper class, SuperMap, which deals with a
  few shortcomings of this class.
*/

#ifndef AMANZI_OPERATORS_SUPER_MAP_HH_
#define AMANZI_OPERATORS_SUPER_MAP_HH_

#include "Teuchos_RCP.hpp"
#include "AmanziComm.hh"

#include "dbc.hh"
#include "Mesh.hh"


namespace Amanzi {

class CompositeVectorSpace;

namespace Operators {

class SuperMapLumped {
 public:
  // Constructor
  SuperMapLumped(const Comm_ptr_type& comm,
           const std::vector<std::string>& compnames,
           const std::vector<int>& dofnums,
           const std::vector<Teuchos::RCP<const Epetra_BlockMap> >& maps,
           const std::vector<Teuchos::RCP<const Epetra_BlockMap> >& ghost_maps);

  SuperMapLumped(const SuperMapLumped& other) = delete;
  virtual ~SuperMapLumped() = default;

  // meta-data
  bool HasComponent(const std::string& compname) const;

  // map accessors
  Teuchos::RCP<const Epetra_Map> Map() const { return map_; }
  Teuchos::RCP<const Epetra_Map> GhostedMap() const { return ghosted_map_; }

  // -- component map accessors
  Teuchos::RCP<const Epetra_BlockMap>
  ComponentMap(const std::string& compname) const {
    return comp_maps_.at(compname);
  }

  Teuchos::RCP<const Epetra_BlockMap>
  ComponentGhostedMap(const std::string& compname) {
    return comp_ghosted_maps_.at(compname);
  }
  
  // index accessors
  const std::vector<int>& Indices(const std::string& compname, int dofnum) const;
  const std::vector<int>& GhostIndices(const std::string& compname, int dofnum) const;

  // block indices.  This is an array of integers, length Map().MyLength(),
  // where each dof and component have a unique integer value.  The returned
  // int is the number of unique values, equal to
  // sum(NumDofs(comp) for comp in components), in this array.
  std::pair<int, Teuchos::RCP<std::vector<int> > > BlockIndices() const;

  // meta-data accessors
  int Offset(const std::string& compname) const { return offsets_.at(compname); }
  int GhostedOffset(const std::string& compname) const { return ghosted_offsets_.at(compname); }
  int NumOwnedElements(const std::string& compname) const { return counts_.at(compname); }
  int NumUsedElements(const std::string& compname) const {
    return counts_.at(compname) + ghosted_counts_.at(compname); }
  int NumDofs(const std::string& compname) const { return num_dofs_.at(compname); }

 protected:
  // iterate over compnames
  typedef std::vector<std::string>::const_iterator name_iterator;
  name_iterator begin() const { return compnames_.begin(); }
  name_iterator end() const { return compnames_.end(); }
  unsigned int size() const { return compnames_.size(); }

 protected:
  // Constructs and returns the vector of indices for a given component into the SuperMapLumped
  virtual const std::vector<int>& CreateIndices_(const std::string& compname, int dofnum, bool ghosted) const;

  // step one of the construction process
  //
  // After this, CreateIndices_() can be called
  void CreateIndexing_();

  // step two of the construction process
  //
  // This creates the SuperMapLumped and uses CreateIndices_() to populate and
  // create the ghosted SuperMapLumped.
  void CreateMap_(const Comm_ptr_type& comm);

 protected:
  std::vector<std::string> compnames_;
  std::map<std::string,int> offsets_;
  std::map<std::string,int> num_dofs_;
  std::map<std::string,int> counts_;
  std::map<std::string,int> ghosted_offsets_;
  std::map<std::string,int> ghosted_counts_;

  int n_local_;
  int n_local_ghosted_;
  
  mutable std::map<std::string, std::map<int, std::vector<int> > > indices_;
  mutable std::map<std::string, std::map<int, std::vector<int> > > ghosted_indices_;

  Teuchos::RCP<Epetra_Map> map_;
  Teuchos::RCP<Epetra_Map> ghosted_map_;

  std::map<std::string, Teuchos::RCP<const Epetra_BlockMap> > comp_maps_;
  std::map<std::string, Teuchos::RCP<const Epetra_BlockMap> > comp_ghosted_maps_;
};


Teuchos::RCP<SuperMapLumped> createSuperMapLumped(const CompositeVectorSpace& cv);

}  // namespace Operators
}  // namespace Amanzi

#endif
