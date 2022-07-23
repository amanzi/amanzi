/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

//! Storage and parameters for BDF1 problems.

/*

  NOTE: Documentation for this file lives in BDF1_TI.hh

*/

#ifndef AMANZI_BDF1STATE_HH_
#define AMANZI_BDF1STATE_HH_

#include <limits>

#include "errors.hh"
#include "dbc.hh"

#include "SolutionHistory.hh"

#include "TimestepController.hh"
#include "TimestepControllerFactory.hh"

namespace Amanzi {

template<class Vector>
struct BDF1_State {
  BDF1_State() :
      freeze_pc(false),
      maxpclag(0),
      extrapolate_guess(true),
      uhist_size(2),
      pc_calls(0),
      seq(-1),
      failed_solve(0),
      failed_current(0),
      pc_updates(0),
      hmin(std::numeric_limits<double>::max()),
      hmax(std::numeric_limits<double>::min()),
      solve_itrs(0) {}

  // Parameters and control
  bool freeze_pc;  // freeze initial preconditioner
  int maxpclag;  // maximum iterations that the preconditioner can be lagged
  bool extrapolate_guess;  // extrapolate forward in time or use previous
                           // step as initial guess for nonlinear solver
  int extrapolation_order;

  // Solution history
  Teuchos::RCP<SolutionHistory<Vector> > uhist;

  // internal counters and state
  int uhist_size;  // extrapolation order for initial guess
  double hlast;  // last step size
  double hpc;  // step size built into the current preconditioner
  int pc_lag;  // counter for how many iterations the preconditioner has been lagged
  int pc_calls;  // counter for the number of preconditioner calls

  // performance counters
  int seq;  // number of steps taken
  int failed_solve;  // number of total failed BCE steps
  int failed_current; // number of current cycle failed steps
  int pc_updates;  // counter for the number of preconditioner updates
  double hmin, hmax;  // minimum and maximum dt used on a successful step

  // performane of nonlinear solver
  int solve_itrs;

  // restart fine constrol
  double tol_multiplier, tol_multiplier_damp;

  void InitializeFromPlist(Teuchos::ParameterList&, const Teuchos::RCP<const Vector>&);
};


/* ******************************************************************
* Initiazition of fundamental parameters
****************************************************************** */
template<class Vector>
void BDF1_State<Vector>::InitializeFromPlist(Teuchos::ParameterList& plist,
        const Teuchos::RCP<const Vector>& initvec) {
  
  // preconditioner control
  freeze_pc = plist.get<bool>("freeze preconditioner", false);
  maxpclag = plist.get<int>("max preconditioner lag iterations", 0);

  // forward time extrapolation (fix me lipnikov@lanl.gov)
  extrapolate_guess = plist.get<bool>("extrapolate initial guess", true);
  extrapolation_order = plist.get<int>("nonlinear iteration initial guess extrapolation order", 1);
  if (extrapolation_order == 0) extrapolate_guess = false;

  // solution history object
  double t0 = plist.get<double>("initial time", 0.0);
  uhist = Teuchos::rcp(new SolutionHistory<Vector>(uhist_size, t0, *initvec));

  // restart fine control
  tol_multiplier = plist.get<double>("restart tolerance relaxation factor", 1.0);
  tol_multiplier_damp = plist.get<double>("restart tolerance relaxation factor damping", 1.0);
}

}  // namespace Amanzi

#endif

