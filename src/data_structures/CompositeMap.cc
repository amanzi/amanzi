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

#include <algorithm>

#include "dbc.hh"
#include "errors.hh"

#include "CompositeMap.hh"

namespace Amanzi {

Import_ptr_type
CompositeMap::Importer(const std::string& name) const
{
  auto import = importers_[name];
  if (import != Teuchos::null) return import;

  auto source_map = ComponentMap(name, false);
  if (source_map == Teuchos::null) return Teuchos::null;

  auto target_map = ComponentMap(name, true);
  if (target_map == Teuchos::null) return Teuchos::null;

  import = Teuchos::rcp(new Import_type(source_map, target_map));
  importers_[name] = import;
  return import;
}



Teuchos::RCP<const CompositeMap>
createCompositeMap(const Comm_ptr_type& comm,
                           const std::vector<std::string>& names,
                           const std::map<std::string, BlockMap_ptr_type>& maps,
                           const std::map<std::string, BlockMap_ptr_type>& ghost_maps,
                           const std::map<std::string, std::size_t> n_dofs)
{
  for (const auto& name : names) {
    if (std::count(names.begin(), names.end(), name) != 1) {
      Errors::Message msg;
      msg << "CompositeMap: duplicate component name: \"" << name << "\"";
      throw(msg);
    }
    
    if (!maps.count(name)) {
      Errors::Message msg;
      msg << "CompositeMap: mis-match between name and map: \"" << name << "\"";
      throw(msg);
    }
    if (!ghost_maps.count(name)) {
      Errors::Message msg;
      msg << "CompositeMap: mis-match between name and ghost_map: \"" << name << "\"";
      throw(msg);
    }
    if (!n_dofs.count(name)) {
      Errors::Message msg;
      msg << "CompositeMap: mis-match between name and number of dofs: \"" << name << "\"";
      throw(msg);
    }
  }
  if (names.size() != maps.size()) {
      Errors::Message msg;
      msg << "CompositeMap: mis-match between name and map.";
      throw(msg);
  }
  if (names.size() != ghost_maps.size()) {
      Errors::Message msg;
      msg << "CompositeMap: mis-match between name and ghost_map.";
      throw(msg);
  }
  if (names.size() != n_dofs.size()) {
      Errors::Message msg;
      msg << "CompositeMap: mis-match between name and number of dofs.";
      throw(msg);
  }
  auto owned_map = Teuchos::rcp(new CompMap(comm, names, maps, n_dofs));
  auto ghosted_map = Teuchos::rcp(new CompMap(comm, names, ghost_maps, n_dofs));
  return Teuchos::rcp(new CompositeMap(owned_map, ghosted_map));
}


Teuchos::RCP<const CompositeMap>
createCompositeMap(const Comm_ptr_type& comm,
                           const std::vector<std::string>& names,
                           const std::vector<BlockMap_ptr_type>& maps,
                           const std::vector<BlockMap_ptr_type>& ghost_maps,
                           const std::vector<std::size_t> n_dofs)
{
  for (const auto& name : names) {
    if (std::count(names.begin(), names.end(), name) != 1) {
      Errors::Message msg;
      msg << "CompositeMap: duplicate component name: \"" << name << "\"";
      throw(msg);
    }
  }
  if (names.size() != maps.size()) {
      Errors::Message msg;
      msg << "CompositeMap: mis-match between name and map.";
      throw(msg);
  }
  if (names.size() != ghost_maps.size()) {
      Errors::Message msg;
      msg << "CompositeMap: mis-match between name and ghost_map.";
      throw(msg);
  }
  if (names.size() != n_dofs.size()) {
      Errors::Message msg;
      msg << "CompositeMap: mis-match between name and number of dofs.";
      throw(msg);
  }
  std::map<std::string, BlockMap_ptr_type> maps_map, ghost_maps_map;
  std::map<std::string, std::size_t> n_dofs_map;
  for (std::size_t i=0; i!=names.size(); ++i) {
    maps_map[names[i]] = maps[i];
    ghost_maps_map[names[i]] = ghost_maps[i];
    n_dofs_map[names[i]] = n_dofs[i];
  }
  auto owned_map = Teuchos::rcp(new CompMap(comm, names, maps_map, n_dofs_map));
  auto ghosted_map = Teuchos::rcp(new CompMap(comm, names, ghost_maps_map, n_dofs_map));
  return Teuchos::rcp(new CompositeMap(owned_map, ghosted_map));
}

} // namespace Amanzi


