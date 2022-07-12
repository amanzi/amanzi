#ifndef MRG_IMIM_TEST_FNBASE_HH_
#define MRG_IMIM_TEST_FNBASE_HH_


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerbosityLevel.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_DataAccess.h"

#include "FnBaseDefs.hh"
#include "ode_2d.hh"
#include <vector>
#include <cmath>
#include "MRG_IMIM_FnBase.hh"


// ODE for testing
class MRG_IMIM_nonlinear_ODE : public Amanzi::ode_2d<Amanzi::MRG_IMIM_FnBase<Epetra_MultiVector>> {
  using Vector = Epetra_MultiVector;
  using Matrix = Epetra_CrsMatrix;


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
  MRG_IMIM_nonlinear_ODE(double atol, double rtol, bool exact_jacobian, 
  Teuchos::RCP< Matrix> A_fast,
  Teuchos::RCP< Matrix> A_slow, Epetra_MpiComm* comm) :
      ode_2d(A_fast, A_slow, comm), exact_jacobian_(exact_jacobian), atol_(atol), rtol_(rtol){

        Epetra_Map map(2,0,*comm);
        temp_memory_ = Teuchos::rcp(new Epetra_Vector(map));
        
        Pu_fast = Teuchos::rcp(new Matrix(*A_));
        Pu_full = Teuchos::rcp(new Matrix(*A_partition_));
        
  }


  

  void FunctionalResidualFull(double t_old, double t_new, std::vector<double> scalings,
  Teuchos::RCP<Vector> u_old_slow,  Teuchos::RCP<Vector> u_new_slow, Teuchos::RCP<Vector> slow_exp,
  Teuchos::RCP<Vector> u_old_fast,  Teuchos::RCP<Vector> u_new_fast, Teuchos::RCP<Vector> fast_exp,
  Teuchos::RCP<Vector> f_slow, Teuchos::RCP<Vector> f_fast) override{
    *f_fast = *u_new_fast;
    *f_slow = *u_new_slow;

    f_fast->Update(1.0, *fast_exp, -1.0);
    f_slow->Update(1.0, *slow_exp, -1.0);

    FunctionalTimeDerivativeFast(t_new, u_new_fast, temp_memory_);

    f_fast->Update(scalings[0], *temp_memory_, 1.0);
    f_slow->Update(scalings[2], *temp_memory_, 1.0);

    FunctionalTimeDerivativeSlow(t_new, u_new_slow, temp_memory_);

    f_fast->Update(scalings[1], *temp_memory_, 1.0);
    f_slow->Update(scalings[3], *temp_memory_, 1.0);

  }

  
  void FunctionalResidualFast(double t_old, double t_new, double scaling,
  Teuchos::RCP<Vector> u_old_fast,  Teuchos::RCP<Vector> u_new_fast,
  Teuchos::RCP<Vector> fast_exp,  Teuchos::RCP<Vector> f_fast){

    *f_fast = *u_new_fast;
    f_fast->Update(1.0, *fast_exp, -1.0);

    FunctionalTimeDerivativeFast(t_new, u_new_fast, temp_memory_);

    f_fast->Update(scaling, *temp_memory_, 1.0);

  }

  // applies preconditioner to u and returns the result in Pu
  int ApplyPreconditionerFull(Teuchos::RCP<const Vector> u_full, Teuchos::RCP<Vector> u_eval) {
    Pu_full->ApplyInverse(*u_full, *u_eval);
  }

  int ApplyPreconditionerFast(Teuchos::RCP<const Vector> u_fast, Teuchos::RCP<Vector> u_eval){
    Pu_fast->ApplyInverse(*u_fast, *u_eval);
  }
  
  // computes a norm on u-du and returns the result
  // Same as standard BE. Any difference needed?
  double ErrorNorm(Teuchos::RCP<const Vector> u, Teuchos::RCP<const Vector> du){
    double norm_du, norm_u;
    du->NormInf(&norm_du);
    u->NormInf(&norm_u);
    return norm_du / ( atol_ + rtol_ * norm_u);
  }

  // updates the preconditioner
  void UpdatePreconditionerFull(double t, std::vector<double> scalings, Teuchos::RCP<const Vector> up_slow, Teuchos::RCP<const Vector> up_fast){
    // do nothing since the preconditioner is the identity

    if (exact_jacobian_) {
      Pu_full -> PutScalar(0.0);
      
      for(auto i = 0; i < matrix_partition_n; ++i) (*Pu_full)[i][i] = -1.0;

      //FIXME: Temporary Setup. Could do something better but not sure how to do block scalings or transformations scalings in Trillinos.
      // A Left scaling and transformation pointer wise to block form would work

      Teuchos::RCP<Matrix> jac_scal = Teuchos::rcp(new Matrix(*Pu_full));
      *jac_scal = *A_partition_;
      
      block_scaling(matrix_n, matrix_n,matrix_n, matrix_n, scalings, jac_scal);

      for (auto i = 0; i < Pu_full->NumMyRows(); i++)
      {
        Pu_full->SumIntoMyValues(i, Pu_full->NumMyCols(), (*jac_scal)[i], indicies);
      }
    } else {
      //TODO: Not implemented yet
    }
  }

  void UpdatePreconditionerFast(double t, double scaling, Teuchos::RCP<const Vector> up){
    
    if (exact_jacobian_) {
      Pu_fast -> PutScalar(0.0);
      
      for(auto i = 0; i < matrix_n; ++i) (*Pu_fast)[i][i] = -1.0;

      for (auto i = 0; i < matrix_n; i++)
      {
        for (auto j = 0; j < matrix_n; j++)
        {
          (*Pu_full)[i][j] += scaling * (*A_fast_)[i][j];
        }
        
      }
      
    } else {
      //TODO: Not implemented yet
    }
  }

  // check the admissibility of a solution
  // override with the actual admissibility check
  // TODO: Not too sure on whether this needs to be done for both fast and slow
   bool IsAdmissible(Teuchos::RCP<const Vector> up){
    return true;
   }

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not

  // TODO: Is full or fast needed for modify predictor (I guess for performance?)
   bool ModifyPredictorFull(double h, Teuchos::RCP<const Vector> u0_fast, Teuchos::RCP<const Vector> u0_slow,
                                   Teuchos::RCP<Vector> u_fast, Teuchos::RCP<Vector> u_slow){
                                    return false;
                                   }

   bool ModifyPredictorFast(double h, Teuchos::RCP<const Vector> u0_fast, Teuchos::RCP<Vector> u_fast) 
   {
    return false;
   }


  // TODO: Is full or fast needed for modify correction (Perhaps for performance?)
   Amanzi::AmanziSolvers::FnBaseDefs::ModifyCorrectionResult 
      ModifyCorrectionFull(double h, Teuchos::RCP<const Vector> res_slow, Teuchos::RCP<const Vector> res_fast,
                       Teuchos::RCP<const Vector> u_slow, Teuchos::RCP<Vector> du_slow,
                       Teuchos::RCP<const Vector> u_fast, Teuchos::RCP<Vector> du_fast) 
   {
    return Amanzi::AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
   }


   Amanzi::AmanziSolvers::FnBaseDefs::ModifyCorrectionResult 
      ModifyCorrectionFast(double h, Teuchos::RCP<const Vector> res,
                       Teuchos::RCP<const Vector> u, Teuchos::RCP<Vector> du)
   {
    return Amanzi::AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
   }


  // update the continuation parameter
   void UpdateContinuationParameter(double lambda) {};

  // calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
   void ChangedSolution(){}

  // experimental routine -- returns the number of linear iterations.
   int ReportStatistics() { return 0; }


  bool exact_jacobian_;
  double atol_, rtol_;
  Teuchos::RCP<Matrix> Pu_fast;
  Teuchos::RCP<Matrix> Pu_full;
   Teuchos::RCP<Vector> temp_memory_;
  const int indicies[4] = {0,1,2,3};
};



#endif
