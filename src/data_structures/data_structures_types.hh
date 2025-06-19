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
enum class InitMode {
  NOALLOC = 0, // construct the data structure, but don't allocate memory
  NONE = 1,    // allocate memory, don't initialize
  COPY = 2,    // allocate, copy from other vector
  ZERO = 3     // allocate, initalize to zero
};

// Amanzi assumes that vector memory is initialized to zero.
//
// Currently Trilinos (Epetra) does this for us, so we don't have to.  Tpetra
// may change that, so we provide the flexibility to change this to ZERO in the
// future.
static const InitMode init_mode_default = InitMode::NONE;

} // namespace Amanzi

#endif
