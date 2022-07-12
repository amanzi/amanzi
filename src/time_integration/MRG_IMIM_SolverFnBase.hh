/*
TODO: Finish up nonlinear Solver Interface


*/

#ifndef AMANZI_MRG_IMIM_SOLVER_FNBASE_
#define AMANZI_MRG_IMIM_SOLVER_FNBASE_

#include "Teuchos_RCP.hpp"

#include "SolverFnBase.hh"
#include "MRG_IMIM_FnBase.hh"
#include "PartitionFnBase.hh"


namespace Amanzi {

template<class Vector>
class MRG_IMIM_SolverFnBase : public AmanziSolvers::SolverFnBase<Vector> {
 public:
 //isfull determines if interface is for full or fast partition
 //TODO: Setup isfull selection
  MRG_IMIM_SolverFnBase(Teuchos::ParameterList& plist,
                    const Teuchos::RCP<MRG_IMIM_FnBase<Vector> >& fn, bool isfull) :
      plist_(plist),
      fn_(fn),
      isfull_(isfull) {};

  // SolverFnBase interface
  // ---------------------------
  // computes the non-linear functional r = F(u)
  virtual void Residual(const Teuchos::RCP<Vector>& u,
                        const Teuchos::RCP<Vector>& r);

  // preconditioner application
  virtual int ApplyPreconditioner(const Teuchos::RCP<const Vector>& r,
                                   const Teuchos::RCP<Vector>& Pr);

  // preconditioner update
  virtual void UpdatePreconditioner(const Teuchos::RCP<const Vector>& u);

  // error norm
  virtual double ErrorNorm(const Teuchos::RCP<const Vector>& u,
                           const Teuchos::RCP<const Vector>& du);

  // Check the admissibility of an inner iterate (ensures preconditions for
  // F(u) to be defined).
  virtual bool IsAdmissible(const Teuchos::RCP<const Vector>& up);

  // Hack a correction for some reason.
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(const Teuchos::RCP<const Vector>& res,
                       const Teuchos::RCP<const Vector>& u,
                       const Teuchos::RCP<Vector>& du);

  // parameter for continuation method
  virtual void UpdateContinuationParameter(double lambda);    

  // bookkeeping for state
  virtual void ChangedSolution();

  //Update the scaling of the Jacobian
  void SetScaling(double scale){
    scaling_ = scale;
  }

  // Other methods
  void SetTimes(double t_old, double t_new) {
    t_old_ = t_old;
    t_new_ = t_new;
    h_ = t_new - t_old;
  }

  void SetPreviousTimeSolution(const Teuchos::RCP<Vector>& u_old) { u_old_ = u_old; }

  void SetExplicitTerms(const Teuchos::RCP<Vector>& u_expterms) {u_expterms_ = u_expterms;}

 protected:
  Teuchos::ParameterList plist_;

  bool isfull_ = 0.0;
  double t_new_ = 0.0;
  double t_old_ = 0.0;
  double h_ = 0.0;
  //TODO: Add options of different scalings: Fast needs to be scalar but full is vector of size 4
  double scaling_ = 0.0;
  Teuchos::RCP<Vector> u_old_;
  Teuchos::RCP<Vector> u_expterms_;

  Teuchos::RCP<MRG_IMIM_FnBase<Vector> > fn_;
};


// computes the non-linear functional r = F(u)
template<class Vector>
void MRG_IMIM_SolverFnBase<Vector>::Residual(const Teuchos::RCP<Vector>& u,
                                         const Teuchos::RCP<Vector>& r) {
  fn_->FastFunctionalResidual(t_old_, t_new_, u_old_, u_expterms_, u, r);
}

// preconditioner application
template<class Vector>
int  MRG_IMIM_SolverFnBase<Vector>::ApplyPreconditioner(const Teuchos::RCP<const Vector>& r,
                                                    const Teuchos::RCP<Vector>& Pr) {
  return fn_->ApplySlowPreconditioner(r, Pr);
}

// preconditioner update
template<class Vector>
void MRG_IMIM_SolverFnBase<Vector>::UpdatePreconditioner(const Teuchos::RCP<const Vector>& u) {
  fn_->UpdateSlowPreconditioner(t_new_, scaling_ * h_, u);
}

// error norm
template<class Vector>
double MRG_IMIM_SolverFnBase<Vector>::ErrorNorm(const Teuchos::RCP<const Vector>& u,
                                            const Teuchos::RCP<const Vector>& du) {
  return fn_->ErrorNorm(u, du);
}


// Check the admissibility of an inner iterate (ensures preconditions for
// F(u) to be defined).
template<class Vector>
bool MRG_IMIM_SolverFnBase<Vector>::IsAdmissible(const Teuchos::RCP<const Vector>& up) {
  return fn_->IsAdmissible(up);
}

// Hack a correction for some reason.
template<class Vector>
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
MRG_IMIM_SolverFnBase<Vector>::ModifyCorrection(const Teuchos::RCP<const Vector>& res,
        const Teuchos::RCP<const Vector>& u,
        const Teuchos::RCP<Vector>& du) {
  return fn_->ModifyCorrection(h_, res, u, du);
}

template<class Vector>
void MRG_IMIM_SolverFnBase<Vector>::UpdateContinuationParameter(double lambda) {
  fn_->UpdateContinuationParameter(lambda);
}

// bookkeeping for state
template<class Vector>
void MRG_IMIM_SolverFnBase<Vector>::ChangedSolution() {
  fn_->ChangedSolution();
}

} // namespace Amanzi
#endif

