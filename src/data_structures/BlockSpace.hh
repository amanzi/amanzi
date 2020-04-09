/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

/*
  Data Structures

  Interface for BlockSpace, a dictionary from keys to their maps.
*/

#ifndef AMANZI_BLOCK_SPACE_HH_
#define AMANZI_BLOCK_SPACE_HH_

#include <vector>
#include <map>
#include "AmanziTypes.hh"
#include "errors.hh"

namespace Amanzi {

class BlockSpace {
 public:
  BlockSpace(const Comm_ptr_type& comm, const std::vector<std::string>& names,
             const std::map<std::string, BlockMap_ptr_type>& master_maps,
             const std::map<std::string, BlockMap_ptr_type>& ghost_maps,
             const std::map<std::string, std::size_t>& num_vectors);
  BlockSpace(const BlockSpace& other) = default;
  BlockSpace& operator=(const BlockSpace&) = default;

  //
  // Map interface
  // ---------------------------------------------
  // equality via approximate SameAs()
  bool SameAs(const BlockSpace& other) const;
  bool LocallySameAs(const BlockSpace& other) const;
  bool SubsetOf(const BlockSpace& other) const;

  Comm_ptr_type Comm() const { return comm_; }

  GO getGlobalLength(bool ghosted = false) const;
  LO getLocalLength(bool ghosted = false) const;

  //
  // Component meta-data
  // ---------------------------------------------
  // Component name information
  bool HasComponent(const std::string& name) const
  {
    return (bool)master_maps_.count(name);
  }

  using name_iterator = std::vector<std::string>::const_iterator;
  name_iterator begin() const { return names_.begin(); }
  name_iterator end() const { return names_.end(); }
  std::size_t size() const { return names_.size(); }

  // accessors
  std::size_t getNumVectors(const std::string& name) const
  {
    return num_vectors_.at(name);
  }
  BlockMap_ptr_type
  ComponentMap(const std::string& name, bool ghosted = false) const;
  Import_ptr_type Importer(const std::string& name) const
  {
    return importers_.at(name);
  }

 protected:
  Comm_ptr_type comm_;
  std::vector<std::string> names_;
  std::map<std::string, BlockMap_ptr_type> master_maps_, ghost_maps_;
  std::map<std::string, std::size_t> num_vectors_;

  std::map<std::string, Import_ptr_type> importers_;
};

} // namespace Amanzi


#endif
