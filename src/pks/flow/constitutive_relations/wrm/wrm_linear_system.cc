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
WRMLinearSystem::WRMLinearSystem(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};

void WRMLinearSystem::InitializeFromPlist_() {
  sat_at_zero_pc_ = plist_.get<double>("saturation at pc=0", 1.0);
  if (plist_.isParameter("alpha")) {
    alpha_ = plist_.get<double>("alpha");
  } else {
    double max_pc = plist_.get<double>("max pc");
    alpha_ = -sat_at_zero_pc_/max_pc;
  }
};

}  // namespace
}  // namespace
}  // namespace
