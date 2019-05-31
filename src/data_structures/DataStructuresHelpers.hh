/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
  ATS

  License: see $ATS_DIR/COPYRIGHT
  Author: Ethan Coon

  Just a few handy typedefs.
------------------------------------------------------------------------- */

#ifndef DATA_STRUCTURE_HELPERS_HH_
#define DATA_STRUCTURE_HELPERS_HH_

#include "AmanziTypes.hh"
#include "AmanziMap.hh"

namespace Amanzi {

// Controls initialization in copy constructor.
typedef enum { INIT_MODE_NONE,
               INIT_MODE_ZERO,
               INIT_MODE_COPY,
               INIT_MODE_NOALLOC} InitMode;

// same as
template<class Map>
bool
SameAs(const Map& one, const Map& two) {
  return one.SameAs(two);
}

template<>
inline bool
SameAs<BlockMap_type>(const BlockMap_type& one, const BlockMap_type& two) {
  return one.locallySameAs(two);
}


} // namespace Amanzi

#endif
