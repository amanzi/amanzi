#ifndef _BDF1STATE_HH_
#define _BDF1STATE_HH_

// This class is based on Neil Carlson's BDF2_DAE module
// that is part of LANL's Truchas code.

#include "dbc.hh"
#include "BDF2_SolutionHistory.hpp"
#include <limits>

namespace Amanzi {

enum bdf_nonlinear_solver_t { BDFNKA, BDFJFNK }; 

struct BDF1State {
  
 public:
  
  BDF1State() {
    seq = -1;
    usable_pc = false;
    mitr = 20;
    ntol = 0.1;

    pcfun_calls = 0;
    updpc_calls = 0;
    updpc_failed = 0;
    retried_bce = 0;
    failed_bce = 0;
    rejected_steps = 0;
    hmin = std::numeric_limits<double>::min();
    hmax = std::numeric_limits<double>::max();
    
    hlimit = 1e10;
    elimit = 1e15;
    maxpclag = 0;

    nonlinear_solver = BDFNKA;

    verbose = false;
  }

  ~BDF1State() {
    if (! uhist) delete uhist;
  }

  void init_solution_history(BDF2::SolutionHistory* sh) {
    ASSERT(sh);
    uhist = sh;
  }


  int       seq;          // number of steps taken
  double    hlast;        // last step size
  double    hpc;          // step size built into the current preconditioner
  bool      usable_pc;    // whether the current preconditioner is usable
  int       mitr;         // maximum number of nonlinear iterations, more and we fail 
  int       minitr;       // minimum number of nonlinear iterations (we will increase time step here)
  int       maxitr;       // maximum number of nonlinear iterations (we cut time step here)
  int       maxpclag;     // maximum iterations that the preconditioner can be lagged
  int       pclagcount;   // counter for how many iterations the preconditioner has been lagged
  double    hlimit;       // maximum allowed time step
  double    elimit;
  double    ntol;         // nonlinear solver error tolerance (relative to 1)
  double    atol, rtol;   // parameters that define the norm to be used in the model evaluator
  BDF2::SolutionHistory* uhist; // solution history structure

  // performance counters
  int    pcfun_calls;    // number of calls to PCFUN
  int    updpc_calls;    // number of calls to UPDPC
  int    updpc_failed;   // number of UPDPC calls returning an error
  int    retried_bce;    // number of retried BCE steps
  int    failed_bce;     // number of completely failed BCE steps
  int    rejected_steps; // number of steps rejected on error tolerance
  double hmin;           // minimum step size used on a successful step
  double hmax;           // maximum step size used on a successful step
  double hinc;           // stepsize increase factor
  double hred;           // stepsize reduction factor

  bdf_nonlinear_solver_t nonlinear_solver;

  // Diagnostics
  bool   verbose;

};

}

#endif // _BDF1STATE_HH_
