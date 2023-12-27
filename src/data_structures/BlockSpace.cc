/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>
/*
  Data Structures

  Implementation for BlockSpace, a dictionary from keys to their maps.
*/

#include "errors.hh"
#include "exceptions.hh"
#include "AmanziMap.hh"
#include "BlockSpace.hh"

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
BlockSpace::getComponentMap(const std::string& name, bool ghosted) const
{
  if (!hasComponent(name)) {
    Errors::Message message("Map: Requested component (" + name + ") does not exist.");
    Exceptions::amanzi_throw(message);
  }
  return ghosted ? ghost_maps_.at(name) : master_maps_.at(name);
}


bool
BlockSpace::locallySameAs(const BlockSpace& other) const
{
  if (names_ != other.names_) return false;
  if (num_vectors_ != other.num_vectors_) return false;
  for (const auto& name : names_) {
    if (!getComponentMap(name)->locallySameAs(*other.getComponentMap(name))) return false;
  }
  return true;
}


bool
BlockSpace::isSameAs(const BlockSpace& other) const
{
  if (names_ != other.names_) return false;
  if (num_vectors_ != other.num_vectors_) return false;
  for (const auto& name : names_) {
    if (!getComponentMap(name)->isSameAs(*other.getComponentMap(name))) return false;
  }
  return true;
}

bool
BlockSpace::isSubsetOf(const BlockSpace& other) const
{
  for (const auto& name : *this) {
    if (!other.hasComponent(name)) return false;
    if (!other.getComponentMap(name, false)->locallySameAs(*getComponentMap(name, false)))
      return false;
  }
  return true;
}


GO
BlockSpace::getGlobalLength(bool ghosted) const
{
  GO count = 0;
  for (const auto& name : names_) {
    count += getComponentMap(name, ghosted)->getGlobalNumElements() * getNumVectors(name);
  }
  return count;
}


LO
BlockSpace::getLocalLength(bool ghosted) const
{
  GO count = 0;
  for (const auto& name : names_) {
    count += getComponentMap(name, ghosted)->getLocalNumElements() * getNumVectors(name);
  }
  return count;
}


} // namespace Amanzi
