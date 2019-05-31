/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

/*
  A dictionary of names to maps.
*/

#ifndef AMANZI_COMPOSITE_MAP_HH_
#define AMANZI_COMPOSITE_MAP_HH_

#include <vector>
#include "Teuchos_RCP.hpp"

#include "dbc.hh"
#include "CompMap.hh"
#include "DataStructuresHelpers.hh"

namespace Amanzi {

template<typename Scalar> class CompositeVector_;

class CompositeMap {

public:
  // Constructor
  CompositeMap(const Teuchos::RCP<const CompMap>& master,
                       const Teuchos::RCP<const CompMap>& ghosted)
      : master_(master), ghosted_(ghosted) {}

  CompositeMap(const CompositeMap& other) = default;
  CompositeMap& operator=(const CompositeMap&) = delete;
  
  // // CompositeMap is a factory of CompositeVectors
  // template<typename Scalar>
  // Teuchos::RCP<CompositeVector_<Scalar> > Create() const {
  //   return Teuchos::rcp(new CompositeVector_<Scalar>(*this));
  // }
    
  // Checks equality
  bool SameAs(const CompositeMap& other, bool modulo_ghosted=false) const {
    if (!Amanzi::SameAs(*master_, *other.master_)) return false;
    if (!modulo_ghosted) return Amanzi::SameAs(*ghosted_, *other.ghosted_);
    return true;
  }

  // CompositeVectors exist on a Comm
  Comm_ptr_type Comm() const { return master_->Comm(); }
  
  // Components are refered to by names.
  // -- Iteration over names of the space
  using name_iterator = CompMap::name_iterator;
  name_iterator begin() const { return master_->begin(); }
  name_iterator end() const { return master_->end(); }
  unsigned int size() const { return master_->size(); }

  std::size_t NumComponents() const {
    return master_->NumComponents();
  }
  bool HasComponent(const std::string& name) const {
    return master_->HasComponent(name);
  }

  // Each component has a number of Degrees of Freedom.
  std::size_t NumVectors(const std::string& name) const {
    return master_->NumVectors(name);
  }

  GO GlobalLength(bool ghosted=false) const {
    return ghosted ? ghosted_->GlobalLength() : master_->GlobalLength();
  }
  
  // ComponentMap access
  BlockMap_ptr_type ComponentMap(const std::string& name, bool ghosted=false) const {
    return ghosted ? ghosted_->ComponentMap(name) : master_->ComponentMap(name);
  }

  // CompMap access
  Teuchos::RCP<const CompMap> Map(bool ghosted=false) const {
    return ghosted ? ghosted_ : master_;
  }

  // Importer access
  Import_ptr_type Importer(const std::string& name) const;
  
 protected:
  Teuchos::RCP<const CompMap> master_;
  Teuchos::RCP<const CompMap> ghosted_;

  mutable std::map<std::string, Import_ptr_type> importers_;
};


Teuchos::RCP<const CompositeMap>
createCompositeMap(const Comm_ptr_type& comm,
                           const std::vector<std::string>& names,
                           const std::vector<BlockMap_ptr_type>& maps,
                           const std::vector<BlockMap_ptr_type>& ghost_maps,
                           const std::vector<std::size_t> n_dofs);

Teuchos::RCP<const CompositeMap>
createCompositeMap(const Comm_ptr_type& comm,
                           const std::vector<std::string>& names,
                           const std::map<std::string, BlockMap_ptr_type>& maps,
                           const std::map<std::string, BlockMap_ptr_type>& ghost_maps,
                           const std::map<std::string, std::size_t> n_dofs);



} // namespace Amanzi

#endif
