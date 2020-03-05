/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_SOLVER_DEFS_HH_
#define AMANZI_SOLVER_DEFS_HH_

namespace Amanzi {
namespace AmanziSolvers {

enum class Monitor {
  UPDATE = 0,
  RESIDUAL = 1,
  //  PRECONDITIONED_RESIDUAL = 2 // NOTE: Never used...  
};

enum class MonitorNorm {
  LINF = 1, // accept decrease in L_INF norm
  L2 = 3,    // accept decrease in L2 norm
  ENORM = 5  // accept decrease in PK's ErrorNorm()
};

enum class MonitorStatus {
  CONTINUE = 1,
  CONVERGED = 0,
  MAX_ITERATIONS = -1,
  DIVERGED = -2,
  STAGNATING = -3,
  DIVERGING = -4,
  INADMISSIBLE_SOLUTION = -5,
  INTERNAL_EXCEPTION = -6,
  BAD_SEARCH_DIRECTION = -7,
  LINEAR_SOLVER_ERROR = -8
};

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
