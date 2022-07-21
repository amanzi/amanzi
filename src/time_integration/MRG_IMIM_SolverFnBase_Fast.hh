/*
TODO: Finish up nonlinear Solver Interface


*/

#ifndef AMANZI_MRG_IMIM_SOLVER_FAST_FNBASE_
#define AMANZI_MRG_IMIM_SOLVER_FAST_FNBASE_

#include "Teuchos_RCP.hpp"

#include "SolverFnBase.hh"
#include "MRG_IMIM_FnBase.hh"
#include "PartitionFnBase.hh"
#include "Epetra_DataAccess.h"


namespace Amanzi {

template<class Vector, class MultiVector>
class MRG_IMIM_SolverFnBase_Fast : public AmanziSolvers::SolverFnBase<Vector> {

 protected:
  Teuchos::ParameterList plist_;
  double t_new_ = 0.0;
  double t_old_ = 0.0;
  double h_ = 0.0;
  double scaling_ = 0.0;
  Teuchos::RCP<Vector> u_old_fast_;
  Teuchos::RCP<Vector> u_expterms_fast_;

  Teuchos::RCP<Amanzi::MRG_IMIM_FnBase<Vector, MultiVector> > fn_;


 public:
 //isfull determines if interface is for full or fast partition
 //TODO: Setup isfull selection
  MRG_IMIM_SolverFnBase_Fast(Teuchos::ParameterList& plist,
                    const Teuchos::RCP<MRG_IMIM_FnBase<Vector, MultiVector> >& fn) :
      plist_(plist),
      fn_(fn) {};

      
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

  void SetPreviousTimeSolution(const Teuchos::RCP<Vector>& u_old_fast) 
  { 
    u_old_fast_ = u_old_fast; 
  }

  void SetExplicitTerms(const Teuchos::RCP<Vector>& u_expterms_fast) {u_expterms_fast_ = u_expterms_fast;}

};

// computes the non-linear functional r = F(u)
template<class Vector, class MultiVector>
void MRG_IMIM_SolverFnBase_Fast<Vector, MultiVector>::Residual(const Teuchos::RCP<Vector>& u,
                                        const Teuchos::RCP<Vector>& r) {
  
  fn_->FunctionalResidualFast(t_old_, t_new_, scaling_, u_old_fast_, u_expterms_fast_, u, r);
  
}

// preconditioner application
template<class Vector, class MultiVector>
int  MRG_IMIM_SolverFnBase_Fast<Vector, MultiVector>::ApplyPreconditioner(const Teuchos::RCP<const Vector>& r,
                                                    const Teuchos::RCP<Vector>& Pr) {
  return fn_->ApplyPreconditionerFast(r, Pr);
}

// preconditioner update
template<class Vector, class MultiVector>
void MRG_IMIM_SolverFnBase_Fast<Vector, MultiVector>::UpdatePreconditioner(const Teuchos::RCP<const Vector>& u) {
  return fn_->UpdatePreconditionerFast(t_new_, scaling_, u);
}

// error norm
template<class Vector, class MultiVector>
double MRG_IMIM_SolverFnBase_Fast<Vector, MultiVector>::ErrorNorm(const Teuchos::RCP<const Vector>& u,
                                            const Teuchos::RCP<const Vector>& du) {
  return fn_->ErrorNorm(u, du);
}


// Check the admissibility of an inner iterate (ensures preconditions for
// F(u) to be defined).
template<class Vector, class MultiVector>
bool MRG_IMIM_SolverFnBase_Fast<Vector, MultiVector>::IsAdmissible(const Teuchos::RCP<const Vector>& up) {
  return fn_->IsAdmissible(up);
}

// Hack a correction for some reason.
template<class Vector, class MultiVector>
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult MRG_IMIM_SolverFnBase_Fast<Vector, MultiVector>::ModifyCorrection(const Teuchos::RCP<const Vector>& res,
                                                                  const Teuchos::RCP<const Vector>& u,
                                                                  const Teuchos::RCP<Vector>& du) {
  return fn_->ModifyCorrectionFast(h_, res, u, du);
}

template<class Vector, class MultiVector>
void MRG_IMIM_SolverFnBase_Fast<Vector, MultiVector>::UpdateContinuationParameter(double lambda) {
  fn_->UpdateContinuationParameter(lambda);
}

// bookkeeping for state
template<class Vector, class MultiVector>
void MRG_IMIM_SolverFnBase_Fast<Vector, MultiVector>::ChangedSolution() {
  fn_->ChangedSolution();
}
} // namespace Amanzi
#endif

