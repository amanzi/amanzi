/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Just a few handy typedefs.
   ------------------------------------------------------------------------- */

#ifndef DATA_STRUCTURE_TYPES_HH_
#define DATA_STRUCTURE_TYPES_HH_

namespace Amanzi {

  // ConstructMode
  // Indicates how copy constructors work.
  //  CONSTRUCT_WITH_NEW_DATA : creates the vector with new uninitialized data
  //  CONSTRUCT_WITH_OLD_DATA : creates a new vector shell with pointers to
  //                            the same old data
  //  CONSTRUCT_WITHOUT_DATA : creates a vector whose CreateData() method has
  //                           not been called
  typedef enum { CONSTRUCT_WITH_NEW_DATA,
                 CONSTRUCT_WITH_OLD_DATA,
                 CONSTRUCT_WITHOUT_DATA } ConstructMode;

} // namespace

#endif
