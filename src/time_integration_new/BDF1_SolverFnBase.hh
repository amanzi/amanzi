/*
  Authors: Ethan Coon (ecoon@lanl.gov)

  Interface that wraps a BDF_FnBase providing the interface for Solver_FnBase
  as used by the BDF1 time integrator.

*/


#ifndef AMANZI_BDF1_SOLVER_FNBASE_
#define AMANZI_BDF1_SOLVER_FNBASE_

#include "Teuchos_RCP.hpp"

#include "SolverFnBase.hh"
#include "BDFFnBase.hh"


namespace Amanzi {

template<class Vec>
class BDF1_SolverFnBase : public AmanziSolvers::SolverFnBase<Vec> {

 public:

  BDF1_SolverFnBase(Teuchos::ParameterList& plist,
                    const Teuchos::RCP<BDFFnBase<Vec> >& bdf_fn) :
      plist_(plist),
      bdf_fn_(bdf_fn) {}

  // SolverFnBase interface
  // ---------------------------
  // computes the non-linear functional r = F(u)
  virtual void Residual(const Teuchos::RCP<Vec>& u,
						const Teuchos::RCP<Vec>& r);

  // preconditioner application
  virtual void ApplyPreconditioner(const Teuchos::RCP<const Vec>& r,
		  const Teuchos::RCP<Vec>& Pr);

  // preconditioner update
  virtual void UpdatePreconditioner(const Teuchos::RCP<const Vec>& u);

  // error norm
  virtual double ErrorNorm(const Teuchos::RCP<const Vec>& u,
                           const Teuchos::RCP<const Vec>& du);

  // Check the admissibility of an inner iterate (ensures preconditions for
  // F(u) to be defined).
  virtual bool IsAdmissible(const Teuchos::RCP<const Vec>& up);

  // Hack a correction for some reason.
  virtual bool ModifyCorrection(const Teuchos::RCP<const Vec>& res,
          const Teuchos::RCP<const Vec>& u,
          const Teuchos::RCP<Vec>& du);

  // bookkeeping for state
  virtual void ChangedSolution();

  // Other methods
  // ---------------------------
  void SetTimes(double t_old, double t_new) {
    t_old_ = t_old;
    t_new_ = t_new;
    h_ = t_new - t_old;
  }

  void SetPreviousTimeSolution(const Teuchos::RCP<Vec>& u_old) { u_old_ = u_old; }

 protected:
  Teuchos::ParameterList plist_;

  double t_new_;
  double t_old_;
  double h_;
  Teuchos::RCP<Vec> u_old_;

  Teuchos::RCP<BDFFnBase<Vec> > bdf_fn_;

};


// SolverFnBase interface
// ---------------------------
// computes the non-linear functional r = F(u)
template<class Vec>
void BDF1_SolverFnBase<Vec>::Residual(const Teuchos::RCP<Vec>& u,
        const Teuchos::RCP<Vec>& r) {
  bdf_fn_->fun(t_old_, t_new_, u_old_, u, r);
}

// preconditioner application
template<class Vec>
void BDF1_SolverFnBase<Vec>::ApplyPreconditioner(const Teuchos::RCP<const Vec>& r,
        const Teuchos::RCP<Vec>& Pr) {
  bdf_fn_->precon(r, Pr);
}

// preconditioner update
template<class Vec>
void BDF1_SolverFnBase<Vec>::UpdatePreconditioner(const Teuchos::RCP<const Vec>& u) {
  bdf_fn_->update_precon(t_new_, u, h_);
}

// error norm
template<class Vec>
double BDF1_SolverFnBase<Vec>::ErrorNorm(const Teuchos::RCP<const Vec>& u,
             const Teuchos::RCP<const Vec>& du) {
  return bdf_fn_->enorm(u, du);
}


// Check the admissibility of an inner iterate (ensures preconditions for
// F(u) to be defined).
template<class Vec>
bool BDF1_SolverFnBase<Vec>::IsAdmissible(const Teuchos::RCP<const Vec>& up) {
  return bdf_fn_->is_admissible(up);
}

// Hack a correction for some reason.
template<class Vec>
bool BDF1_SolverFnBase<Vec>::ModifyCorrection(const Teuchos::RCP<const Vec>& res,
        const Teuchos::RCP<const Vec>& u,
        const Teuchos::RCP<Vec>& du) {
  return bdf_fn_->modify_correction(h_, res, u, du);
}

// bookkeeping for state
template<class Vec>
void BDF1_SolverFnBase<Vec>::ChangedSolution() {
  bdf_fn_->changed_solution();
}


} // namespace
#endif
