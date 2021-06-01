/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT

MPC for the 1D lake model.

------------------------------------------------------------------------- */
#include "mpc_lake_1D.hh"

namespace Amanzi {

RegisteredPKFactory<MPCLake1D> MPCLake1D::reg_("coupled lake-soil model");

} // namespace
