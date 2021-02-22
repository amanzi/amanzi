/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Basic land cover/plant function type

#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {


LandCover::LandCover(Teuchos::ParameterList& plist) :
  rooting_depth_max(plist.get<double>("rooting depth max [m]")),
  rooting_profile_alpha(plist.get<double>("rooting profile alpha [-]")),
  rooting_profile_beta(plist.get<double>("rooting profile beta [-]")),
  stomata_closed_mafic_potential(plist.get<double>("mafic potential at fully closed stomata [Pa]")),
  stomata_open_mafic_potential(plist.get<double>("mafic potential at fully open stomata [Pa]")),
  mannings_n(plist.get<double>("Manning's n [?]")),
  leaf_on_doy(plist.get<double>("leaf on time [doy]")),
  leaf_off_doy(plist.get<double>("leaf off time [doy]")),
  albedo_ground(plist.get<double>("albedo of ground surface [-]")),
  emissivity_ground(plist.get<double>("emissivity of ground surface [-]")),
  albedo_canopy(plist.get<double>("albedo of canopy [-]")),
  emissivity_canopy(plist.get<double>("emissivity of canopy [-]")),
  beers_k_sw(plist.get<double>("Beer's law extinction coefficient, shortwave [-]")),
  beers_k_lw(plist.get<double>("Beer's law extinction coefficient, longwave [-]")),
  snow_transition_depth(plist.get<double>("snow transition depth [m]")),
  dessicated_zone_thickness(plist.get<double>("dessicated zone thickness [m]")),
  clapp_horn_b(plist.get<double>("Clapp and Hornberger b [-]"))
{}


LandCoverMap getLandCover(Teuchos::ParameterList& plist)
{
  LandCoverMap lc;
  for (auto& item : plist) {
    if (plist.isSublist(item.first)) {
      lc.insert({item.first, LandCover{plist.sublist(item.first)}});
    }
  }
  return lc;
}



} // namespace SurfaceBalance
} // namespace ATS
