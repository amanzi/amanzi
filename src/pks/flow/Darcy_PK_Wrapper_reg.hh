// Darcy_PK_Wrapper registration
#include "Darcy_PK_Wrapper.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Darcy_PK_Wrapper> Darcy_PK_Wrapper::reg_("Darcy");

} // namespace
} // namespace
