// Richards_PK_Wrapper registration
#include "Richards_PK_Wrapper.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Richards_PK_Wrapper> Richards_PK_Wrapper::reg_("richards");

} // namespace
} // namespace
