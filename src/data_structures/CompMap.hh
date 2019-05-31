/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Data Structures

  Interface for CompMap, a dictionary from keys to their maps.
*/

#ifndef AMANZI_COMP_MAP_HH_
#define AMANZI_COMP_MAP_HH_

#include <vector>
#include <map>
#include "AmanziTypes.hh"
#include "errors.hh"

namespace Amanzi {

class CompMap {
 public:
  CompMap(const Comm_ptr_type& comm,
          const std::vector<std::string> names,
          const std::map<std::string, BlockMap_ptr_type>& maps,
          const std::map<std::string, std::size_t>& num_dofs)
      : comm_(comm),
        names_(names),
        maps_(maps),
        num_dofs_(num_dofs) {}

  CompMap(const CompMap& other) = delete;
  CompMap& operator=(const CompMap&) = delete;

  // equality via approximate SameAs()
  bool SameAs(const CompMap& other) const;
  
  // iterate over names
  typedef std::vector<std::string>::const_iterator name_iterator;
  name_iterator begin() const { return names_.begin(); }
  name_iterator end() const { return names_.end(); }
  unsigned int size() const { return names_.size(); }

  Comm_ptr_type Comm() const { return comm_; }
  
  // sizes
  GO GlobalLength() const;
  std::size_t NumComponents() const { return size(); }
  std::size_t NumVectors(const std::string& name) const { return num_dofs_.at(name); }

  // accessors
  bool HasComponent(const std::string& name) const { return (bool) maps_.count(name); }

  BlockMap_ptr_type ComponentMap(const std::string& name) const {
    if (!HasComponent(name)) {
      Errors::Message message("Map: Requested component ("+name+") does not exist.");
      throw(message);
    }
    return maps_.at(name);
  }

 protected:
  Comm_ptr_type comm_;
  std::vector<std::string> names_;
  std::map<std::string, BlockMap_ptr_type> maps_;
  std::map<std::string, std::size_t> num_dofs_;
};

} // namespace Amanzi


#endif
