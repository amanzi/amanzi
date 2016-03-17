
#ifndef AMANZI_SOLVER_NOX_
#define AMANZI_SOLVER_NOX_

#include "Teuchos_RCP.hpp"

#include "NOX_Solver_Factory.H"

#include "ResidualDebugger.hh"
#include "SolverFnBase.hh"
#include "NoxSolverFnBase.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector, class VectorSpace>
class SolverNox : public Solver<Vector,VectorSpace> {
 public:
  SolverNox(Teuchos::ParameterList& plist) :
      plist_(plist) {}

  // still needs tests!  
  virtual void Init(const Teuchos::RCP<SolverFnBase<Vector> >& fn,
                    const VectorSpace& map) {
    group_ = Teuchos::rcp(new NoxSolverFnBase<Vector>(fn));
    solver_ = NOX::Solver::buildSolver(group, tests, plist_);
  }

  virtual int Solve(const Teuchos::RCP<Vector>& u) {
    NoxVector u_nox(u);
    group_->setX(u_nox);
    solver_->solve()
  }

  // mutators
  virtual void set_tolerance(double tol) { // needs to pass to status test
  }

  virtual void set_pc_lag(double pc_lag) {}
  virtual void set_db(const Teuchos::RCP<ResidualDebugger>& db) {}
  
  // access 
  virtual double tolerance() = 0;
  virtual double residual() = 0;
  virtual int num_itrs() { solver_->getNumIterations(); }
  virtual int returned_code() = 0;
  virtual int pc_calls() = 0;
  virtual int pc_updates() = 0;


 private:
  Teuchos::RCP<NoxSolverFnBase<Vector> > group_;
  Teuchos::RCP<NOX::Solver::Generic> solver_;
  Teuchos::ParameterList plist_;
};

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
