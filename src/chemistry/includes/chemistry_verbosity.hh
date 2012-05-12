/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_VERBOSITY_HH_
#define AMANZI_CHEMISTRY_VERBOSITY_HH_

/*!
**  \file chemistry_verbosity.hh
**
**  \author Ben Andre
**  \date 2011-09
**
**  General description of the verbosity functionality is:
**
**  VerbosityFlags is an an array that controls whether a particular
**  verbosity level is enabled.
**
**  Verbosity is an enumeration of indicies into the VerbosityFlags array
**
**  VerbosityMap is a map that converts a string (from an input file),
**  into a Verbosity enum index value.
**
**  Use:
**
**  The user supplies a string for a verbosity level. Convert it into
**  an index, then "set" the corresponding value in the flags
**  array. Later, check if that flag is set using the "test".
**
**  std::string verbosity_level = strings::kDebugDatabase;
**  
**  VerbosityFlags verbosity_flags;
**  VerbosityMap verbosity_map = CreateVerbosityMap();
**  int index = verbosity_map.at(verbosity_level);
**  verbosity_flags.set(index);
**  if (verbosity_flags.test(kDebugDatabase)) { do something }
**
**  TODO(bandre): need some logic that if the silent flag is set,then
**  everything else is ignored....?
**
*/
#include <iostream>
#include <string>
#include <map>
#include <bitset>

#include "chemistry_strings.hh"

namespace amanzi {
namespace chemistry {

enum Verbosity { 
  kSilent,
  kTerse,
  kVerbose,
  kError,
  kWarning,
  kDebug,
  kDebugDriver,
  kDebugInputFile,
  kDebugDatabase,
  kDebugActivityModel,
  kDebugSpeciation,
  kDebugLinearSolver,
  kDebugChemistryProcessKernel,
  kDebugBeaker,
  kDebugMineralKinetics,
  // old stuff is indented
                 kDebugSorptionIsotherm,
                 kDebugIonExchange,
                 kDebugNever  // always last!
};

typedef std::map<std::string, Verbosity> VerbosityMap;

typedef std::bitset<32> VerbosityFlags;

VerbosityMap CreateVerbosityMap(void);

}  // namespace chemistry
}  // namespace amanzi
#endif     /* AMANZI_CHEMISTRY_VERBOSITY_HH_ */
