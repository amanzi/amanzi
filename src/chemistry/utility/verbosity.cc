/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*******************************************************************************
 **
 **  Helper functions for verbosity
 **
 *******************************************************************************/
#include "verbosity.hh"

#include <string>
#include <map>

namespace amanzi {
namespace chemistry {

VerbosityMap CreateVerbosityMap(void)
{
  // create a std::map to convert a string into a verbosity enum value

  VerbosityMap verbosity_map;
  verbosity_map[kSilentName] = kSilent;
  verbosity_map[kTerseName] = kTerse;
  verbosity_map[kVerboseName] = kVerbose;
  verbosity_map[kDebugName] = kDebug;
  verbosity_map[kDebugBeakerName] = kDebugBeaker;
  return verbosity_map;
}  // end CreateVerbosityMap()

}  // namespace chemistry
}  // namespace amanzi
