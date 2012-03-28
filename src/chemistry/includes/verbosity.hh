/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_VERBOSITY_HH_
#define AMANZI_CHEMISTRY_VERBOSITY_HH_

#include <string>
#include <map>

namespace amanzi {
namespace chemistry {

enum Verbosity { kSilent,
                 kTerse,
                 kVerbose,
                 kDebug,
                 kDebugBeaker,
                 kDebugInputFile,
                 kDebugMineralKinetics,
                 kDebugSorptionIsotherm,
                 kDebugIonExchange,
                 kDebugNewtonSolver,
                 kDebugChemistryProcessKernel,
                 kDebugNever  // always last!
};

static const std::string kSilentName("silent");
static const std::string kTerseName("terse");
static const std::string kVerboseName("verbose");
static const std::string kDebugName("debug");
static const std::string kDebugBeakerName("debug_beaker");
static const std::string kDebugInputFileName("debug_input_file");

typedef std::map<std::string, Verbosity> VerbosityMap;

VerbosityMap CreateVerbosityMap(void);

}  // namespace chemistry
}  // namespace amanzi
#endif     /* AMANZI_CHEMISTRY_VERBOSITY_HH_ */
