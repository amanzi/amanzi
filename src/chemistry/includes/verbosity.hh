/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_VERBOSITY_HH_
#define AMANZI_CHEMISTRY_VERBOSITY_HH_

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


}  // namespace chemistry
}  // namespace amanzi
#endif     /* AMANZI_CHEMISTRY_VERBOSITY_HH_ */
