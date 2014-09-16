/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the derived MPCStrong class.  Is both a PK and a Model
Evalulator, providing needed methods for BDF time integration of the coupled
system.

Completely automated and generic to any sub PKs, this uses a block diagonal
preconditioner.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */


#include "MPCStrong.hh"

namespace Amanzi {

template<>
RegisteredPKFactory<MPCStrong<ImplicitFnPK> > MPCStrong<ImplicitFnPK>::reg_("strong MPC");

} // namespace
