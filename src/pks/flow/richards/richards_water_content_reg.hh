#include "richards_water_content.hh"
namespace Amanzi {
namespace Flow {
Utils::RegisteredFactory<FieldEvaluator,RichardsWaterContent> RichardsWaterContent::reg_("richards water content");
}
}
