#ifndef MRG_EXIM_TEST_FNBASE_HH_
#define MRG_EXIM_TEST_FNBASE_HH_


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerbosityLevel.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_DataAccess.h"

#include "MRG_EXIM_FnBase.hh"
#include "FnBaseDefs.hh"
#include "ode_2d.hh"
#include <vector>
#include <cmath>

// ODE for testing
class MRG_EXIM_linear2d_ODE :  public Amanzi::ode_2d<Amanzi::MRG_EXIM_FnBase<Epetra_MultiVector>>  {
  using Vector = Epetra_MultiVector;
  using Matrix = Epetra_CrsMatrix;

private:

  bool exact_jacobian_;
  double atol_, rtol_;
  Teuchos::RCP<Matrix> Pu_slow;
  Teuchos::RCP<Vector> temp_memory_;


public:

  /**
   * @brief Construct a new mgii nonlinear ODE object
   * 
   * Problem is defined as follows
   * 
   * y' = A y
   * 
   * where y \in R^{2} and A \in R^{2 x 2}
   * 
   * README: Does this need a Eptra_MPIComm ? It is just a test.
   * 
   * @param atol 
   * @param rtol 
   * @param exact_jacobian 
   * @param A 
   */
  MRG_EXIM_linear2d_ODE(double atol, double rtol, bool exact_jacobian, 
  Teuchos::RCP< Matrix> A_fast, Teuchos::RCP< Matrix> A_slow, Epetra_MpiComm* comm) :
      ode_2d(A_fast, A_slow, comm), exact_jacobian_(exact_jacobian), atol_(atol), rtol_(rtol){

        Epetra_Map map(2,0,*comm);
        temp_memory_ = Teuchos::rcp(new Epetra_Vector(map));
        
        Pu_slow = Teuchos::rcp(new Matrix(*A_));
        
  }

  void SlowFunctionalResidual(double t_old, double t_new, double scaling,
  Teuchos::RCP<Vector> u_old,  Teuchos::RCP<Vector> u_exp,
  Teuchos::RCP<Vector> u_new,  Teuchos::RCP<Vector> f_eval) override {

    *f_eval = *u_new;
    f_eval->Update(1.0, *u_old , 1.0, *u_exp, -1.0);


    FunctionalTimeDerivativeSlow(t_new, u_new, temp_memory_);

    f_eval->Update(1.0, *temp_memory_, 1.0/scaling);

  }

  int  ApplySlowPreconditioner(Teuchos::RCP<const Vector> u_slow, Teuchos::RCP<Vector> u_eval) override {
    return Pu_slow->ApplyInverse(*u_slow, *u_eval);
  }
  
  // computes a norm on u-du and returns the result
  double ErrorNorm(Teuchos::RCP<const Vector> u, Teuchos::RCP<const Vector> du) override {
    double norm_du, norm_u;
    du->NormInf(&norm_du);
    u->NormInf(&norm_u);
    return norm_du / ( atol_ + rtol_ * norm_u);
  }

  /**
   * @brief Produces the Preconditioner for solving the nonlinear system
   * 
   *  The preconditioner is:   I/scaling + J(t, u)
   * 
   * @param t 
   * @param scaling 
   * @param up 
   */
  void UpdateSlowPreconditioner(double t, double scaling, Teuchos::RCP<const Vector> u) override{
    
    if (exact_jacobian_) {
      Pu_slow -> PutScalar(0.0);

      for (auto i = 0; i < matrix_n; ++i)
      {
        for (auto j = 0; j < matrix_n; ++j)
        {
          (*Pu_slow)[i][j] = (*A_slow_)[i][j];

          if (i == j) (*Pu_slow)[i][j] -= 1.0 / scaling;
        }
        
      }
      
    } else {
      //TODO: Not implemented yet
    }
  }

  // check the admissibility of a solution
  // override with the actual admissibility check
  // TODO: Not too sure on whether this needs to be done for both fast and slow
  bool IsAdmissible(Teuchos::RCP<const Vector> up) override {
    return true;
   }

  bool ModifyPredictorSlow(double h, Teuchos::RCP<const Vector> u0_slow, Teuchos::RCP<Vector> u_slow) override
   {
    return false;
   }


  Amanzi::AmanziSolvers::FnBaseDefs::ModifyCorrectionResult 
      ModifyCorrectionSlow(double h, Teuchos::RCP<const Vector> res,
                       Teuchos::RCP<const Vector> u, Teuchos::RCP<Vector> du) override
   {
    return Amanzi::AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
   };


  // update the continuation parameter
   void UpdateContinuationParameter(double lambda) {};

  // calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  void ChangedSolution() override {};

  // experimental routine -- returns the number of linear iterations.
  int ReportStatistics() override { return 0; };

};



#endif
