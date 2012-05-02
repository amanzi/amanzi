/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Markus Berndt (berndt@lanl.gov) 
  Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "wrm_linear_system.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<WRM,WRMLinearSystem> WRMLinearSystem::factory_("linear system");

/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
WRMLinearSystem::WRMLinearSystem(Teuchos::ParameterList& wrm_plist) :
    wrm_plist_(wrm_plist) {
  InitializeFromPlist_();
};

void WRMLinearSystem::InitializeFromPlist_() {
  sat_at_zero_pc_ = wrm_plist_.get<double>("saturation at pc=0", 1.0);
  if (wrm_plist_.isParameter("alpha")) {
    alpha_ = wrm_plist_.get<double>("alpha");
  } else {
    double max_pc = wrm_plist_.get<double>("max pc");
    alpha_ = -sat_at_zero_pc_/max_pc;
  }
};

}  // namespace
}  // namespace
}  // namespace
