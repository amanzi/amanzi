/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Data Structures

  Implementation for BlockSpace, a dictionary from keys to their maps.
*/

#include "BlockSpace.hh"
#include "AmanziMap.hh"

namespace Amanzi {

BlockSpace::BlockSpace(const Comm_ptr_type& comm,
                 const std::vector<std::string>& names,
                 const std::map<std::string, BlockMap_ptr_type>& master_maps,
                 const std::map<std::string, BlockMap_ptr_type>& ghost_maps,
                 const std::map<std::string, std::size_t>& num_vectors)
  : comm_(comm),
    names_(names),
    master_maps_(master_maps),
    ghost_maps_(ghost_maps),
    num_vectors_(num_vectors)
{
  for (const auto& name : names_) {
    if (ghost_maps_.at(name) != Teuchos::null) {
      importers_[name] = Teuchos::rcp(new Import_type(master_maps_.at(name), ghost_maps_.at(name)));
    }
  }
}

BlockMap_ptr_type
BlockSpace::ComponentMap(const std::string& name, bool ghosted) const
{
  if (!HasComponent(name)) {
    Errors::Message message("Map: Requested component ("+name+") does not exist.");
    throw(message);
  }
  return ghosted ? ghost_maps_.at(name) : master_maps_.at(name);
}


bool BlockSpace::LocallySameAs(const BlockSpace& other) const
{
  if (names_ != other.names_) return false;
  if (num_vectors_ != other.num_vectors_) return false;
  for (const auto& name : names_) {
    if (!ComponentMap(name)->locallySameAs(*other.ComponentMap(name))) return false;
  }
  return true;
}


bool BlockSpace::SameAs(const BlockSpace& other) const
{
  if (names_ != other.names_) return false;
  if (num_vectors_ != other.num_vectors_) return false;
  for (const auto& name : names_) {
    if (!ComponentMap(name)->isSameAs(*other.ComponentMap(name))) return false;
  }
  return true;
}


GO BlockSpace::GlobalLength(bool ghosted) const
{
  GO count = 0;
  for (const auto& name : names_) {
    count += ComponentMap(name, ghosted)->getGlobalNumElements() * NumVectors(name);
  }
  return count;
}


LO BlockSpace::MyLength(bool ghosted) const
{
  GO count = 0;
  for (const auto& name : names_) {
    count += ComponentMap(name, ghosted)->getNodeNumElements() * NumVectors(name);
  }
  return count;
}



} // namespace Amanzi
