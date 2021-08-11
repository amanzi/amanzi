/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT

MPC for the 1D lake model (lake+soil+richards).

------------------------------------------------------------------------- */
#include "mpc_lake_soil_richards.hh"

namespace Amanzi {

RegisteredPKFactory<MPCLakeSoilRichards> MPCLakeSoilRichards::reg_("coupled lake-soil-richards model");

} // namespace
