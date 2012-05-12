/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*******************************************************************************
 **
 **  Helper functions for verbosity
 **
 *******************************************************************************/
#include "chemistry_verbosity.hh"

#include <string>
#include <map>

#include "chemistry_utilities.hh"
#include "chemistry_strings.hh"

namespace amanzi {
namespace chemistry {

VerbosityMap CreateVerbosityMap(void)
{
  // create a std::map to convert a string into a verbosity enum value
  // add mixed case and all lower case keys.

  VerbosityMap verbosity_map;
  std::string key;

  verbosity_map[strings::kVerbositySilent] = kSilent;
  utilities::LowerCaseString(strings::kVerbositySilent, &key);
  verbosity_map[key] = kSilent;

  verbosity_map[strings::kVerbosityTerse] = kTerse;
  utilities::LowerCaseString(strings::kVerbosityTerse, &key);
  verbosity_map[key] = kTerse;

  verbosity_map[strings::kVerbosityVerbose] = kVerbose;
  utilities::LowerCaseString(strings::kVerbosityVerbose, &key);
  verbosity_map[key] = kVerbose;

  verbosity_map[strings::kVerbosityWarning] = kWarning;
  utilities::LowerCaseString(strings::kVerbosityWarning, &key);
  verbosity_map[key] = kWarning;

  verbosity_map[strings::kVerbosityError] = kError;
  utilities::LowerCaseString(strings::kVerbosityError, &key);
  verbosity_map[key] = kError;

  verbosity_map[strings::kVerbosityDebugChemistryCoordinator] = kDebugChemistryCoordinator;
  utilities::LowerCaseString(strings::kVerbosityDebugChemistryCoordinator, &key);
  verbosity_map[key] = kDebugChemistryCoordinator;

  verbosity_map[strings::kVerbosityDebugInputFile] = kDebugInputFile;
  utilities::LowerCaseString(strings::kVerbosityDebugInputFile, &key);
  verbosity_map[key] = kDebugInputFile;

  verbosity_map[strings::kVerbosityDebugDatabase] = kDebugDatabase;
  utilities::LowerCaseString(strings::kVerbosityDebugDatabase, &key);
  verbosity_map[key] = kDebugDatabase;

  verbosity_map[strings::kVerbosityDebugReactionNetwork] = kDebugReactionNetwork;
  utilities::LowerCaseString(strings::kVerbosityDebugReactionNetwork, &key);
  verbosity_map[key] = kDebugReactionNetwork;

  verbosity_map[strings::kVerbosityDebugActivityModel] = kDebugActivityModel;
  utilities::LowerCaseString(strings::kVerbosityDebugActivityModel, &key);
  verbosity_map[key] = kDebugActivityModel;

  verbosity_map[strings::kVerbosityDebugSpeciation] = kDebugSpeciation;
  utilities::LowerCaseString(strings::kVerbosityDebugSpeciation, &key);
  verbosity_map[key] = kDebugSpeciation;

  verbosity_map[strings::kVerbosityDebugOperatorSplitting] = kDebugOperatorSplitting;
  utilities::LowerCaseString(strings::kVerbosityDebugOperatorSplitting, &key);
  verbosity_map[key] = kDebugOperatorSplitting;

  verbosity_map[strings::kVerbosityDebugNonlinearSolver] = kDebugNonlinearSolver;
  utilities::LowerCaseString(strings::kVerbosityDebugNonlinearSolver, &key);
  verbosity_map[key] = kDebugNonlinearSolver;

  verbosity_map[strings::kVerbosityDebugODESolver] = kDebugODESolver;
  utilities::LowerCaseString(strings::kVerbosityDebugODESolver, &key);
  verbosity_map[key] = kDebugODESolver;

  verbosity_map[strings::kVerbosityDebugLinearSolver] = kDebugLinearSolver;
  utilities::LowerCaseString(strings::kVerbosityDebugLinearSolver, &key);
  verbosity_map[key] = kDebugLinearSolver;

  verbosity_map[strings::kVerbosityDebugChemistryProcessKernel] = kDebugChemistryProcessKernel;
  utilities::LowerCaseString(strings::kVerbosityDebugChemistryProcessKernel, &key);
  verbosity_map[key] = kDebugChemistryProcessKernel;

  return verbosity_map;
}  // end CreateVerbosityMap()

}  // namespace chemistry
}  // namespace amanzi
