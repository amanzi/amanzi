/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
//! Basic land cover/plant function type

#include "exceptions.hh"
#include "errors.hh"

#include "LandCover.hh"
#include "seb_nan.hh"

namespace Amanzi {
namespace SurfaceBalance {


LandCover::LandCover(Teuchos::ParameterList& plist) :
  rooting_depth_max(plist.get<double>("rooting depth max [m]", NAN)),
  rooting_profile_alpha(plist.get<double>("rooting profile alpha [-]", NAN)),
  rooting_profile_beta(plist.get<double>("rooting profile beta [-]", NAN)),
  stomata_closed_mafic_potential(plist.get<double>("mafic potential at fully closed stomata [Pa]", NAN)),
  stomata_open_mafic_potential(plist.get<double>("mafic potential at fully open stomata [Pa]", NAN)),
  mannings_n(plist.get<double>("Manning's n [?]", NAN)),
  leaf_on_doy(plist.get<double>("leaf on time [doy]", NAN)),
  leaf_off_doy(plist.get<double>("leaf off time [doy]", NAN)),
  albedo_ground(plist.get<double>("albedo of ground surface [-]", NAN)),
  emissivity_ground(plist.get<double>("emissivity of ground surface [-]", NAN)),
  albedo_canopy(plist.get<double>("albedo of canopy [-]", NAN)),
  emissivity_canopy(plist.get<double>("emissivity of canopy [-]", NAN)),
  beers_k_sw(plist.get<double>("Beer's law extinction coefficient, shortwave [-]", NAN)),
  beers_k_lw(plist.get<double>("Beer's law extinction coefficient, longwave [-]", NAN)),
  snow_transition_depth(plist.get<double>("snow transition depth [m]", NAN)),
  water_transition_depth(plist.get<double>("water transition depth [m]", NAN)),
  dessicated_zone_thickness(plist.get<double>("dessicated zone thickness [m]", NAN)),
  clapp_horn_b(plist.get<double>("Clapp and Hornberger b [-]", NAN)),
  bare_ground_surface_roughness(plist.get<double>("roughness length of bare ground [m]", NAN)),
  snow_surface_roughness(plist.get<double>("roughness length of snow-covered ground [m]", NAN))
{}


LandCoverMap getLandCover(Teuchos::ParameterList& plist)
{
  LandCoverMap lc;
  for (auto& item : plist) {
    if (plist.isSublist(item.first)) {
      lc.insert({item.first, LandCover{plist.sublist(item.first)}});
    }
  }
  if (lc.size() == 0) {
    Errors::Message message("LandCover is used, but no entries were found in the 'state->initial conditions->land cover types' list.");
    Exceptions::amanzi_throw(message);
  }
  return lc;
}



} // namespace SurfaceBalance
} // namespace ATS
