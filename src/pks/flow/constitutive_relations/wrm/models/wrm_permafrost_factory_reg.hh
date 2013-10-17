/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for WRM implementations.
   ------------------------------------------------------------------------- */

#include "wrm_permafrost_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// explicity instantitate the static data of Factory<WRM>
//template<> Factory<WRM>::map_type* Factory<WRM>::map_;
template<> Utils::Factory<WRMPermafrostModel>::map_type* Utils::Factory<WRMPermafrostModel>::map_;

} // namespace
} // namespace
} // namespace

