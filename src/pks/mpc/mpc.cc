/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Implementation for the Base MPC class.  A multi process coordinator is a PK
   (process kernel) which coordinates several PKs.  Each of these coordinated PKs
   may be MPCs themselves, or physical PKs.  Note this does NOT provide a full
   implementation of PK -- it does not supply the advance() method.  Therefore
   this class cannot be instantiated, but must be inherited by derived classes
   which finish supplying the functionality.  Instead, this provides the data
   structures and methods (which may be overridden by derived classes) for
   managing multiple PKs.

   Most of these methods simply loop through the coordinated PKs, calling their
   respective methods.
   ------------------------------------------------------------------------- */

#include "mpc.hh"

namespace Amanzi {


} // namespace
