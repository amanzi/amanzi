/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
  ATS

  Just a few handy typedefs.
------------------------------------------------------------------------- */

#ifndef DATA_STRUCTURE_TYPES_HH_
#define DATA_STRUCTURE_TYPES_HH_

namespace Amanzi {

// Controls initialization in copy constructor.
typedef enum { INIT_MODE_NONE, INIT_MODE_ZERO, INIT_MODE_COPY, INIT_MODE_NOALLOC } InitMode;

} // namespace Amanzi

#endif
