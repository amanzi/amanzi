/*
TODO: Finish up nonlinear Solver Interface


*/

#ifndef AMANZI_MRG_IMIM_SOLVER_FULL_FNBASE_
#define AMANZI_MRG_IMIM_SOLVER_FULL_FNBASE_

#include "Teuchos_RCP.hpp"
#include <vector>

#include "SolverFnBase.hh"
#include "MRG_IMIM_FnBase.hh"
#include "PartitionFnBase.hh"
#include "Epetra_DataAccess.h"


namespace Amanzi {

template<class Vector, class MultiVector>
class MRG_IMIM_SolverFnBase_Full : public AmanziSolvers::SolverFnBase<MultiVector> {

    
 protected:
  Teuchos::ParameterList plist_;
  double t_new_ = 0.0;
  double t_old_ = 0.0;
  double h_ = 0.0;
  std::vector<double> scaling_vec_ = 0.0;
  Teuchos::RCP<MultiVector> u_old_full_;
  Teuchos::RCP<MultiVector> u_expterms_full_;

  Teuchos::RCP<Amanzi::MRG_IMIM_FnBase<Vector, MultiVector> > fn_;

 public:
 //isfull determines if interface is for full or fast partition
 //TODO: Setup isfull selection
  MRG_IMIM_SolverFnBase_Full(Teuchos::ParameterList& plist,
                    const Teuchos::RCP<MRG_IMIM_FnBase<Vector, MultiVector> >& fn) :
      plist_(plist),
      fn_(fn),
      scaling_vec_(4) {};

  void SetScaling(std::vector<double> scale_vec){
    scaling_vec_ = scale_vec;
  }

  // Other methods
  void SetTimes(double t_old, double t_new) {
    t_old_ = t_old;
    t_new_ = t_new;
    h_ = t_new - t_old;
  }

  void SetPreviousTimeSolution(const Teuchos::RCP<MultiVector>& u_old_full) 
  { 
    u_old_full_ = u_old_full; 
  }

  void SetExplicitTerms(const Teuchos::RCP<MultiVector>& u_expterms_full) {u_expterms_full_ = u_expterms_full;}

    // SolverFnBase interface
  // ---------------------------
  // computes the non-linear functional r = F(u)
  virtual void Residual(const Teuchos::RCP<MultiVector>& u,
                        const Teuchos::RCP<MultiVector>& r);

  // preconditioner application
  virtual int ApplyPreconditioner(const Teuchos::RCP<const MultiVector>& r,
                                   const Teuchos::RCP<MultiVector>& Pr);

  // preconditioner update
  virtual void UpdatePreconditioner(const Teuchos::RCP<const MultiVector>& u);

  // error norm
  virtual double ErrorNorm(const Teuchos::RCP<const MultiVector>& u,
                           const Teuchos::RCP<const MultiVector>& du);

  // Check the admissibility of an inner iterate (ensures preconditions for
  // F(u) to be defined).
  virtual bool IsAdmissible(const Teuchos::RCP<const MultiVector>& up);

  // Hack a correction for some reason.
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(const Teuchos::RCP<const MultiVector>& res,
                       const Teuchos::RCP<const MultiVector>& u,
                       const Teuchos::RCP<MultiVector>& du);

  // parameter for continuation method
  virtual void UpdateContinuationParameter(double lambda);    

  // bookkeeping for state
  virtual void ChangedSolution();

  


};

// computes the non-linear functional r = F(u)
template<class Vector, class MultiVector>
void MRG_IMIM_SolverFnBase_Full<Vector, MultiVector>::Residual(const Teuchos::RCP<MultiVector>& u,
                                        const Teuchos::RCP<MultiVector>& r){

fn_->FunctionalResidualFull(t_old_, t_new_, scaling_vec_, u_old_full_, u_expterms_full_, u, r);

}

// preconditioner application
template<class Vector, class MultiVector>
int MRG_IMIM_SolverFnBase_Full<Vector, MultiVector>::ApplyPreconditioner(const Teuchos::RCP<const MultiVector>& r,
                                                    const Teuchos::RCP<MultiVector>& Pr) {
return fn_->ApplyPreconditionerFull(r, Pr);
}

// preconditioner update
template<class Vector, class MultiVector>
void MRG_IMIM_SolverFnBase_Full<Vector, MultiVector>::UpdatePreconditioner(const Teuchos::RCP<const MultiVector>& u) {
return fn_->UpdatePreconditionerFull(t_new_, scaling_vec_, u);
}

// error norm
template<class Vector, class MultiVector>
double MRG_IMIM_SolverFnBase_Full<Vector, MultiVector>::ErrorNorm(const Teuchos::RCP<const MultiVector>& u,
                                            const Teuchos::RCP<const MultiVector>& du) {
return fn_->ErrorNorm(u, du);
}


// Check the admissibility of an inner iterate (ensures preconditions for
// F(u) to be defined).
template<class Vector, class MultiVector>
bool MRG_IMIM_SolverFnBase_Full<Vector, MultiVector>::IsAdmissible(const Teuchos::RCP<const MultiVector>& up) {
return fn_->IsAdmissible(up);
}

// Hack a correction for some reason.
template<class Vector, class MultiVector>
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult MRG_IMIM_SolverFnBase_Full<Vector, MultiVector>::ModifyCorrection(const Teuchos::RCP<const MultiVector>& res,
                                                                const Teuchos::RCP<const MultiVector>& u,
                                                                const Teuchos::RCP<MultiVector>& du) {
return fn_->ModifyCorrectionFull(h_, res, u, du);
}

template<class Vector, class MultiVector>
void MRG_IMIM_SolverFnBase_Full<Vector, MultiVector>::UpdateContinuationParameter(double lambda) {
fn_->UpdateContinuationParameter(lambda);
}

// bookkeeping for state
template<class Vector, class MultiVector>
void MRG_IMIM_SolverFnBase_Full<Vector, MultiVector>::ChangedSolution() {
fn_->ChangedSolution();
}



} // namespace Amanzi
#endif

