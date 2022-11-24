/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
  ATS

  License: see $ATS_DIR/COPYRIGHT
  Author: Ethan Coon

  Just a few handy typedefs.
------------------------------------------------------------------------- */

#ifndef DATA_STRUCTURE_TYPES_HH_
#define DATA_STRUCTURE_TYPES_HH_

namespace Amanzi {

// Controls initialization in copy constructor.
typedef enum { INIT_MODE_NONE, INIT_MODE_ZERO, INIT_MODE_COPY, INIT_MODE_NOALLOC } InitMode;

} // namespace Amanzi

#endif
