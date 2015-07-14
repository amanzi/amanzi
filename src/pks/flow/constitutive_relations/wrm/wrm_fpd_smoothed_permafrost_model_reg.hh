/*
Author: Ethan Coon

Painter's permafrost model with freezing point depression.

 */

#include "wrm.hh"
#include "wrm_fpd_smoothed_permafrost_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {


// registry of method
Utils::RegisteredFactory<WRMPermafrostModel,WRMFPDSmoothedPermafrostModel> WRMFPDSmoothedPermafrostModel::factory_("fpd smoothed permafrost model");

} // namespace FlowRelations
} // namespace Flow
} // namespace Amanzi
