/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Basic land cover/plant function type

#include "LandCover.hh"

namespace ATS {
namespace SurfaceBalance {


LandCover::LandCover(Teuchos::ParameterList& plist) :
  max_rooting_depth(plist.get<double>("max rooting depth [m]", -1)),
  rooting_profile_alpha(plist.get<double>("rooting profile alpha [-]", -1)),
  rooting_profile_beta(plist.get<double>("rooting profile beta [-]", -1)),
  mannings_n(plist.get<double>("Manning's n [-]", -1)),
  leaf_on_doy(plist.get<double>("leaf on time [doy]", -1)),
  leaf_off_doy(plist.get<double>("leaf off time [doy]", -1)),
  interception_coef(plist.get<double>("interception coefficient [-]", -1)) {}


LandCoverMap getLandCover(Teuchos::ParameterList& plist)
{
  LandCoverMap lc;
  for (auto& item : plist) {
    if (lc.isSublist(item.first)) {
      lc.try_insert(item.first, LandCover{plist.sublist(item.first)});
    }
  }
  return lc;
}



} // namespace SurfaceBalance
} // namespace ATS
