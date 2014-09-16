#ifndef AMANZI_FNBASE_DEFS_HH_
#define AMANZI_FNBASE_DEFS_HH_

namespace Amanzi {
namespace AmanziSolvers {

// contained in its own namespace for easier using

namespace FnBaseDefs {
// enum for ModifyCorrection control
enum ModifyCorrectionResult {
  CORRECTION_NOT_MODIFIED = 0,
  CORRECTION_MODIFIED = 1,
  CORRECTION_MODIFIED_LAG_BACKTRACKING
};
}

}
}

#endif
