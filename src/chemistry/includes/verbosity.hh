/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_VERBOSITY_HH_
#define AMANZI_CHEMISTRY_VERBOSITY_HH_

enum Verbosity { kSilent,
                 kTerse,
                 kVerbose,
                 kDebug,
                 kDebugBeaker,
                 kDebugInputFile,
                 kDebugMineralKinetics,
                 kDebugIonExchange,
                 kDebugNewtonSolver,
                 kDebugChemistryProcessKernel,
                 kDebugNever  // always last!
};


#endif     /* AMANZI_CHEMISTRY_VERBOSITY_HH_ */
