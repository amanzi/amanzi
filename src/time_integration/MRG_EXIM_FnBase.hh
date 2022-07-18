#ifndef AMANZI_MRG_EXIM_FNBASE_HH_
#define AMANZI_MRG_EXIM_FNBASE_HH_

#include "FnBaseDefs.hh"
#include "PartitionFnBase.hh"

namespace Amanzi {

// This is the interface definition for the BDF class.
// The nonlinear functional, preconditioner, and error
// functions must be derived from this class to be 
// usable with BDF1_TI.

template<class Vector>
class MRG_EXIM_FnBase : public PartitionFnBase<Vector> {
 public:
  virtual ~MRG_EXIM_FnBase() = default;

    /**
   * @brief Evaluate the nonlinear-Functional F(u) = f_eval
   * 
   * The nonlinear functional is 
   * F(u) = 1/scaling(u_exp + u_old - u_new) + f(u_new)
   * 
   * @param t_old 
   * @param t_new 
   * @param scaling dt * a_{i,i}
   * @param u_old 
   * @param u_exp 
   * @param u_new 
   * @param f_eval output
   */
  virtual void SlowFunctionalResidual(double t_old, double t_new, double scaling,
  Teuchos::RCP<Vector> u_old,  Teuchos::RCP<Vector> u_exp,
  const Teuchos::RCP<Vector> u_new,  const Teuchos::RCP<Vector> &f_eval) = 0  ;

  virtual int ApplySlowPreconditioner(const Teuchos::RCP<const Vector>& u_slow, const Teuchos::RCP<Vector>& u_eval) = 0;


  // updates the preconditioner
  virtual void UpdateSlowPreconditioner(double t, double scaling, Teuchos::RCP<const Vector> u) = 0;

  // check the admissibility of a solution
  // override with the actual admissibility check
  virtual bool IsAdmissible(Teuchos::RCP<const Vector> up) = 0;

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not
  virtual bool ModifyPredictorSlow(double h, Teuchos::RCP<const Vector> u0_slow, Teuchos::RCP<Vector> u_slow) = 0;

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult 
      ModifyCorrectionSlow(double h, Teuchos::RCP<const Vector> res,
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
