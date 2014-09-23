/* -------------------------------------------------------------------------
Amanzi

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the derived MPCWeak class.  Provides only the advance()
method missing from MPC.hh.  In weak coupling, we simply loop over the
sub-PKs, calling their advance() methods and returning failure if any fail.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#include "MPCWeak.hh"

namespace Amanzi {

RegisteredPKFactory<MPCWeak> MPCWeak::reg_("weak MPC");

} // namespace
