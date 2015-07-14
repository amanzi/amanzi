/*
Author: Ethan Coon

Painter's permafrost model.

 */

#include "wrm.hh"
#include "wrm_old_permafrost_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

Utils::RegisteredFactory<WRMPermafrostModel,WRMOldPermafrostModel> WRMOldPermafrostModel::factory_("old permafrost model");

} // namespace
} // namespace
} // namespace
