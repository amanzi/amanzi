/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! OutputUtils: Utility functions for I/O.

#include "AmanziTypes.hh"
#include "AmanziVector.hh"

namespace Amanzi {

//
// Given a ghosted, unordered GIDs map, get the ghosted, naturally ordered map.
Map_ptr_type
GetNaturalMap(const Map_ptr_type& ghosted_map, const Map_ptr_type& owned_map)
{
  // create the owned natural map.
  Map_type owned_natural(owned_map->getGlobalNumElements(),
                         owned_map->getLocalNumElements(),
                         0,
                         owned_map->getComm());

  Vector_type_<Map_type::global_ordinal_type> natural(ghosted_map);
  {
    auto nv = natural.getLocalViewHost(Tpetra::Access::ReadWrite);
    for (int i = 0; i != owned_map->getLocalNumElements(); ++i) {
      nv(i, 0) = owned_natural.getGlobalElement(i);
    }
  }
  Import_type importer(owned_map, ghosted_map);
  natural.doImport(natural, importer, Tpetra::INSERT);

  Teuchos::Array<Map_type::global_ordinal_type> inds(natural.getLocalLength());
  natural.get1dCopy(inds);
  auto ghosted_natural = Teuchos::rcp(new Map_type(
    ghosted_map->getGlobalNumElements(), inds, 0, owned_map->getComm()));
  return ghosted_natural;
}

} // namespace Amanzi
