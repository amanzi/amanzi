/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ben Andre

  Verbosity keywords.
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

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
