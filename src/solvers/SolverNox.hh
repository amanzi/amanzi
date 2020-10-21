/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

//! Calls Nox nonlinear solvers/JFNK.


  
#ifndef AMANZI_SOLVER_NOX_
#define AMANZI_SOLVER_NOX_

#include "AmanziGroup.hh"
#include "Teuchos_RCP.hpp"

#include "NOX_Solver_Factory.H"
#include "NOX_Solver_Generic.H"
#include "NOX_StatusTest_Factory.H"
#include "NOX_Utils.H"

#include "MatrixJF.hh"
#include "Matrix.hh"
#include "ResidualDebugger.hh"
#include "SolverFnBase.hh"
#include "Solver.hh"
#include "SolverFnBaseJF.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class VectorClass, class VectorSpace>
class SolverNox : public Solver<VectorClass, VectorSpace> {
 public:
  SolverNox(Teuchos::ParameterList& plist)
  {
    plist_ = Teuchos::rcp(new Teuchos::ParameterList(plist));
    tol_ = 1.e-6;

    // Grab the nonlinear solver portion of the parameterList
    Teuchos::ParameterList& nonlinearList = plist_->sublist("nonlinear solver");

    // Create the parameterList for NOX
    noxList_ = Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& printing = noxList_->sublist("Printing");
    printing.set("Output Information",511);

    // Get the parameters from that
    tol_ = nonlinearList.get("nonlinear tolerance", tol_);
    double overflow_tol = nonlinearList.get("diverged tolerance", 1.0e10);
    int max_itrs = nonlinearList.get("limit iterations", 20);
    int max_divergence_count = nonlinearList.get("max divergent iterations", 3);

    // Create the printing utility
    NOX::Utils utils(NOX::Utils::Debug);

    // Create the convergence tests
    Teuchos::ParameterList stl;
    stl.set("Test Type", "Combo");
    stl.set("Combo Type", "OR");
    stl.set("Number of Tests", 3);
    Teuchos::ParameterList& conv = stl.sublist("Test 0");
    Teuchos::ParameterList& divergence = stl.sublist("Test 1");
    Teuchos::ParameterList& maxiters = stl.sublist("Test 2");
    conv.set("Test Type", "NormF");
    conv.set("Tolerance", tol_);
    conv.set("Norm Type", "Max Norm");
    conv.set("Scale Type", "Unscaled");
    divergence.set("Test Type", "Divergence");
    divergence.set("Tolerance", overflow_tol);
    divergence.set("Consecutive Iterations", max_divergence_count);
    maxiters.set("Test Type", "MaxIters");
    maxiters.set("Maximum Iterations", max_itrs);
    tests_ = NOX::StatusTest::buildStatusTests(stl,utils);
  }

  // still needs tests!
  void Init(const Teuchos::RCP<SolverFnBase<VectorClass> >& fn,
            const VectorSpace& map) {
    // wrap the base fn in a JF fn
    jf_fnbase_ = Teuchos::rcp(new SolverFnBaseJF<VectorClass,VectorSpace>(*plist_, fn, map));

    group_ = Teuchos::rcp(new AmanziGroup<VectorClass>(jf_fnbase_));
  }

  int Solve(const Teuchos::RCP<VectorClass>& u) {
    NoxVector<VectorClass> u_nox(u);
    group_->setX(u_nox);
    // We can't build the solver until we have set X
    solver_ = NOX::Solver::buildSolver(group_, tests_, noxList_);
    solver_->solve();
    *u = *dynamic_cast<const NoxVector<VectorClass>&>(group_->getX()).get_vector();
    return 0;
  }

  // mutators
  void set_tolerance(double tol) { tol_ = tol; }
  void set_pc_lag(double pc_lag) { pc_lag_ = pc_lag; }
  void set_db(const Teuchos::RCP<ResidualDebugger>& db) { db_ = db; }
  
  // access 
  double tolerance() { return tol_; }
  double residual() { return 0.0; }  // FIXME
  int num_itrs() { return solver_->getNumIterations(); }
  int returned_code() { return 0; }
  int pc_calls() { return 0; }
  int pc_updates() { return 0; }

 private:
  Teuchos::RCP<AmanziGroup<VectorClass> > group_;
  Teuchos::RCP<NOX::Solver::Generic> solver_;
  Teuchos::RCP<NOX::StatusTest::Generic> tests_;
  Teuchos::RCP<Teuchos::ParameterList> plist_;
  Teuchos::RCP<Teuchos::ParameterList> noxList_;
  Teuchos::RCP<ResidualDebugger> db_;
  Teuchos::RCP<SolverFnBaseJF<VectorClass,VectorSpace> > jf_fnbase_;
  double tol_;
  double pc_lag_;
};

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
