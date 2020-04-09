/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#ifndef DATA_STRUCTURE_HELPERS_HH_
#define DATA_STRUCTURE_HELPERS_HH_

#include "AmanziTypes.hh"
#include "AmanziMap.hh"

namespace Amanzi {

// Controls initialization in copy constructor.
enum class InitMode { NONE, ZERO, COPY, NOALLOC };

// same as
template <class Map>
bool
SameAs(const Map& one, const Map& two)
{
  return one.SameAs(two);
}

template <>
inline bool
SameAs<BlockMap_type>(const BlockMap_type& one, const BlockMap_type& two)
{
  return one.isSameAs(two);
}


// same as
template <class Map>
bool
LocallySameAs(const Map& one, const Map& two)
{
  return one.LocallySameAs(two);
}

template <>
inline bool
LocallySameAs<BlockMap_type>(const BlockMap_type& one, const BlockMap_type& two)
{
  return one.locallySameAs(two);
}

} // namespace Amanzi

#endif
