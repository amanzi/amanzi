/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef VERBOSITY_HH_

#define VERBOSITY_HH_

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


#endif     /* VERBOSITY_HH_ */

