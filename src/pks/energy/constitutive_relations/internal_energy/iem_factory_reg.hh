/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for IEM implementations.
   ------------------------------------------------------------------------- */

#include "iem_factory.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

// explicity instantitate the static data of Factory<IEM>
//template<> Factory<IEM>::map_type* Factory<IEM>::map_;
template<> Utils::Factory<IEM>::map_type* Utils::Factory<IEM>::map_;

} // namespace
} // namespace
} // namespace

