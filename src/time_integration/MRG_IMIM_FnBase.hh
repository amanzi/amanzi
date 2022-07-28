#ifndef AMANZI_MRG_IMIM_FNBASE_HH_
#define AMANZI_MRG_IMIM_FNBASE_HH_

#include "FnBaseDefs.hh"
#include "PartitionFnBase.hh"
#include <vector>

namespace Amanzi {

/**
 * @brief Interface for the nonlinear solvers on Multi-Rate methods
 * 
 * This provides an interface for PKs and MPKs to provide residual operators for
 * the multi rate time integrators
 * 
 * This interface distincts them into full and fast interfaces
 * 
 * TODO: See if this is wanted to be generalized to support FIRK or others potientially
 * This would be generalized by making the input vectors of vectors and providing coefficents for scaling
 * This also simplifies the inteface
 * 
 * TODO: Potientally could look into using Vector instead of multiple parameters.
 * This would generalize the setup. Also would generlize the nonlinear solver setup.
 * Will have to use full vector due to how the nonlinear solvers are setup
 * 
 * @tparam Vector templated class to provide based of TreeVector.hh
 */

template<class Vector>
class MRG_IMIM_FnBase : public PartitionFnBase<Vector>  {

 public:
  virtual ~MRG_IMIM_FnBase() = default;
  
  // computes the non-linear functional f = f(t,u,udot)
  virtual void FunctionalResidualFull(double t_old, double t_new, std::vector<double> scalings,
      Teuchos::RCP<Vector> u_old_full, Teuchos::RCP<Vector> u_exp_full, 
       const Teuchos::RCP<Vector> u_new_full, const Teuchos::RCP<Vector> &f_eval_full) = 0;

  
  virtual void FunctionalResidualFast(double t_old, double t_new, double scaling,
    Teuchos::RCP<Vector> u_old_fast,  Teuchos::RCP<Vector> u_exp_fast,
    const Teuchos::RCP<Vector> u_new_fast, const Teuchos::RCP<Vector> &f_eval_fast) = 0;

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditionerFull(Teuchos::RCP<const Vector> u_full, Teuchos::RCP<Vector> u_eval) = 0;

  virtual int ApplyPreconditionerFast(Teuchos::RCP<const Vector> u_fast, Teuchos::RCP<Vector> u_eval) = 0;
  
  // updates the preconditioner
  // TODO: Provide interface for coefficents to scale for both Full and Fast
  virtual void UpdatePreconditionerFull(double t, std::vector<double> scalings, Teuchos::RCP<const Vector> u_full) = 0;

  virtual void UpdatePreconditionerFast(double t, double scaling, Teuchos::RCP<const Vector> u_fast) = 0;

  
  // computes a norm on u-du and returns the result
  virtual double ErrorNorm(Teuchos::RCP<const Vector> u, Teuchos::RCP<const Vector> du) = 0;


  // // computes a norm on u-du and returns the result
  // virtual double ErrorNorm(Teuchos::RCP<const Vector> u, Teuchos::RCP<const Vector> du) = 0;

  // check the admissibility of a solution
  // override with the actual admissibility check
  // TODO: Not too sure on whether this needs to be done for both fast and slow
  virtual bool IsAdmissible(Teuchos::RCP<const Vector> u) = 0;

  // check the admissibility of a solution
  // override with the actual admissibility check
  // TODO: Not too sure on whether this needs to be done for both fast and slow
  // virtual bool IsAdmissible(Teuchos::RCP<const Vector> u) = 0;

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not

  // TODO: Is full or fast needed for modify predictor (I guess for performance?)
  virtual bool ModifyPredictorFull(double h, Teuchos::RCP<const Vector> u0_full, Teuchos::RCP<Vector> u_full) = 0;

  virtual bool ModifyPredictorFast(double h, Teuchos::RCP<const Vector> u0_fast, Teuchos::RCP<Vector> u_fast) = 0;


  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated

  // TODO: IS full or fast needed for modify correction (Perhaps for performance?)
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult 
      ModifyCorrectionFull(double h, Teuchos::RCP<const Vector> res_full,
                       Teuchos::RCP<const Vector> u_full, Teuchos::RCP<Vector> du_full) = 0;

  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult 
      ModifyCorrectionFast(double h, Teuchos::RCP<const Vector> res,
                       Teuchos::RCP<const Vector> u, Teuchos::RCP<Vector> du) = 0;

  // update the continuation parameter
  virtual void UpdateContinuationParameter(double lambda) {};

  // calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() = 0;

  // experimental routine -- returns the number of linear iterations.
  virtual int ReportStatistics() { return 0; }
};

}  // namespace Amanzi

#endif 
