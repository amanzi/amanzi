/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ben Andre

  General description of the verbosity functionality is:

  VerbosityFlags is an an array that controls whether a particular
  verbosity level is enabled.
  Verbosity is an enumeration of indicies into the VerbosityFlags array
*/

#ifndef AMANZI_CHEMISTRY_VERBOSITY_HH_
#define AMANZI_CHEMISTRY_VERBOSITY_HH_

#include <iostream>
#include <string>
#include <map>
#include <bitset>

#include "chemistry_strings.hh"

namespace Amanzi {
namespace AmanziChemistry {

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

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
