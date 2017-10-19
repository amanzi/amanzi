#include "mpc_predictor_corrector.hh"

namespace Amanzi {

template<>
RegisteredPKFactory<MPCPredictorCorrector<PK_BDF_Default> > MPCPredictorCorrector<PK_BDF_Default>::reg_("predictor-corrector MPC");

} // namespace
