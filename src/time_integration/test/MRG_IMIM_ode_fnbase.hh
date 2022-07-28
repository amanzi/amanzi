#ifndef MRG_IMIM_TEST_FNBASE_HH_
#define MRG_IMIM_TEST_FNBASE_HH_


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerbosityLevel.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_DataAccess.h"
#include "Epetra_LinearProblem.h"

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
#include "Amesos2_MatrixAdapter.hpp"
#include "Amesos2_SolverCore.hpp"
#include "Amesos2_MultiVecAdapter.hpp"

#include "FnBaseDefs.hh"
#include "ode_2d_Tree.hh"
#include <vector>
#include <cmath>
#include "MRG_IMIM_FnBase.hh"


// ODE for testing
class MRG_IMIM_linear2d_ODE : public Amanzi::ode_2d_Tree<Amanzi::MRG_IMIM_FnBase<Amanzi::TreeVector>> {
  using Vector = Amanzi::TreeVector;
  using Matrix = Epetra_CrsMatrix;

private:

  bool exact_jacobian_;
  double atol_, rtol_;
  Teuchos::RCP<Matrix> Pu_fast;
  Teuchos::RCP<Matrix> Pu_full;
  Teuchos::RCP<Vector> temp_memory_;
  Teuchos::RCP<Epetra_MultiVector> temp_memory_full_x_;
  Teuchos::RCP<Epetra_MultiVector> temp_memory_full_b_;

  //Viewing vectors
  Teuchos::RCP<Vector> view_u_vec_0;
  Teuchos::RCP<Vector> view_u_vec_1;
  Teuchos::RCP<Vector> view_f_vec_0;
  Teuchos::RCP<Vector> view_f_vec_1;

  Teuchos::RCP<Amesos2::Solver<Matrix, Epetra_MultiVector>> solver_fast;
  Teuchos::RCP<Amesos2::Solver<Matrix, Epetra_MultiVector>> solver_coupled;


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
  MRG_IMIM_linear2d_ODE(double atol, double rtol, bool exact_jacobian, 
  Teuchos::RCP< Matrix> A_fast,
  Teuchos::RCP< Matrix> A_slow, Epetra_MpiComm* comm, Teuchos::RCP<Vector> init_vector) :
      ode_2d_Tree(A_fast, A_slow, comm), exact_jacobian_(exact_jacobian), atol_(atol), rtol_(rtol){

        Epetra_Map map(2,0,*comm);
        Teuchos::RCP<Epetra_Vector> init_solver_memory = Teuchos::rcp(new Epetra_Vector(map));

        temp_memory_ = Teuchos::rcp(new Vector(*init_vector));

        Epetra_Map map_full(4,0,*comm);
        temp_memory_full_x_ = Teuchos::rcp(new Epetra_MultiVector(map_full, 1));
        temp_memory_full_b_ = Teuchos::rcp(new Epetra_MultiVector(map_full, 1));
        
        Pu_fast = Teuchos::rcp(new Matrix(*A_));
        Pu_full = Teuchos::rcp(new Matrix(*A_partition_));

        solver_fast = Amesos2::create<Matrix, Epetra_MultiVector>("KLU2", Pu_fast, init_solver_memory, init_solver_memory);
        solver_coupled = Amesos2::create<Matrix, Epetra_MultiVector>("KLU2", Pu_full);
        
        
  }


  
  void FunctionalResidualFull(double t_old, double t_new, std::vector<double> scalings,
      Teuchos::RCP<Vector> u_old_full, Teuchos::RCP<Vector> u_exp_full, 
       const Teuchos::RCP<Vector> u_new_full, const Teuchos::RCP<Vector> &f_eval_full) override{
    *f_eval_full = *u_new_full;

    f_eval_full->Update(1.0, *u_old_full, 1.0, *u_exp_full, -1.0);
    
    //Apply scalings
    FunctionalTimeDerivativeFast(t_new, u_new_full->SubVector(0), temp_memory_);

    f_eval_full->SubVector(0)->Update(scalings[0], *temp_memory_, 1.0);
    f_eval_full->SubVector(1)->Update(scalings[2], *temp_memory_, 1.0);

    FunctionalTimeDerivativeSlow(t_new, u_new_full->SubVector(1), temp_memory_);

    f_eval_full->SubVector(0)->Update(scalings[1], *temp_memory_, 1.0);
    f_eval_full->SubVector(1)->Update(scalings[3], *temp_memory_, 1.0);
  }

  
  void FunctionalResidualFast(double t_old, double t_new, double scaling,
    Teuchos::RCP<Vector> u_old_fast,  Teuchos::RCP<Vector> u_exp_fast,
    const Teuchos::RCP<Vector> u_new_fast, const Teuchos::RCP<Vector> &f_eval_fast) override{

    *f_eval_fast = *u_new_fast;
    f_eval_fast->Update(1.0, *u_old_fast , 1.0, *u_exp_fast, -1.0);

    FunctionalTimeDerivativeFast(t_new, u_new_fast, temp_memory_);

    f_eval_fast->Update(scaling, *temp_memory_, 1.0);

  }

  // applies preconditioner to u and returns the result in Pu
  int ApplyPreconditionerFull(Teuchos::RCP<const Vector> u_full, Teuchos::RCP<Vector> u_eval) override{
    //FIXME: This just does copying vectors over. In reality this should be reference based using the operators
    //temporary setup

    (*temp_memory_full_b_)[0][0] = (*u_full->SubVector(0)->Data()->ViewComponent("cell"))[0][0];
    (*temp_memory_full_b_)[0][1] = (*u_full->SubVector(0)->Data()->ViewComponent("cell"))[0][1];
    (*temp_memory_full_b_)[0][2] = (*u_full->SubVector(1)->Data()->ViewComponent("cell"))[0][0];
    (*temp_memory_full_b_)[0][3] = (*u_full->SubVector(1)->Data()->ViewComponent("cell"))[0][1];

    solver_coupled->setB(temp_memory_full_b_);
    solver_coupled->setX(temp_memory_full_x_);

    //Cannot do multivectors block solve will get memory leak
    solver_coupled->symbolicFactorization().numericFactorization().solve();

    (*u_eval->SubVector(0)->Data()->ViewComponent("cell"))[0][0] = (*temp_memory_full_x_)[0][0];
    (*u_eval->SubVector(0)->Data()->ViewComponent("cell"))[0][1] = (*temp_memory_full_x_)[0][1];
    (*u_eval->SubVector(1)->Data()->ViewComponent("cell"))[0][0] = (*temp_memory_full_x_)[0][2];
    (*u_eval->SubVector(1)->Data()->ViewComponent("cell"))[0][1] = (*temp_memory_full_x_)[0][3];

    return 0;
  }

  int ApplyPreconditionerFast(Teuchos::RCP<const Vector> u_fast, Teuchos::RCP<Vector> u_eval) override{
    solver_fast->setB(u_fast->Data()->ViewComponent("cell"));
    solver_fast->setX(u_eval->Data()->ViewComponent("cell"));

    solver_fast->symbolicFactorization().numericFactorization().solve();

    return 0;
  }
  
  // computes a norm on u-du and returns the result
  // Same as standard BE. Any difference needed?
  double ErrorNorm(Teuchos::RCP<const Vector> u, Teuchos::RCP<const Vector> du) override{
    double norm_du, norm_u;
    du->NormInf(&norm_du);
    u->NormInf(&norm_u);
    return norm_du;
  }  

  // updates the preconditioner
  void UpdatePreconditionerFull(double t, std::vector<double> scalings, Teuchos::RCP<const Vector> u_full) override{

    if (exact_jacobian_) {
      //FIXME: Temporary Setup. Could do something better but not sure how to do block scalings or transformations scalings in Trillinos.

      *Pu_full = *A_partition_;
      
      block_scaling(matrix_n, matrix_n,matrix_n, matrix_n, scalings, Pu_full);

      for (auto i = 0; i < matrix_partition_n; i++) (*Pu_full)[i][i] -= 1;
      
    } else {
      //TODO: Not implemented yet
    }
  }

  void UpdatePreconditionerFast(double t, double scaling, Teuchos::RCP<const Vector> u_fast) override{
    
    if (exact_jacobian_) {

      for (auto i = 0; i < matrix_n; ++i)
      {
        for (auto j = 0; j < matrix_n; ++j)
        {
          (*Pu_fast)[i][j] = scaling * (*A_fast_)[i][j];

          if (i == j) (*Pu_fast)[i][j] -= 1.0;
        }
        
      }
      
    } else {
      //TODO: Not implemented yet
    }
  }

  // check the admissibility of a solution
  // override with the actual admissibility check
  // TODO: Not too sure on whether this needs to be done for both fast and slow
   bool IsAdmissible(Teuchos::RCP<const Vector> u) override{
    return true;
   }

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not

  // TODO: Is full or fast needed for modify predictor (I guess for performance?)
   bool ModifyPredictorFull(double h, Teuchos::RCP<const Vector> u0_full, Teuchos::RCP<Vector> u_full) override{
                                    return false;
                                   }

   bool ModifyPredictorFast(double h, Teuchos::RCP<const Vector> u0_fast, Teuchos::RCP<Vector> u_fast) override{
    return false;
   }


  // TODO: Is full or fast needed for modify correction (Perhaps for performance?)
   Amanzi::AmanziSolvers::FnBaseDefs::ModifyCorrectionResult 
      ModifyCorrectionFull(double h, Teuchos::RCP<const Vector> res_full,
                       Teuchos::RCP<const Vector> u_full, Teuchos::RCP<Vector> du_full) override
   {
    return Amanzi::AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
   }


   Amanzi::AmanziSolvers::FnBaseDefs::ModifyCorrectionResult 
      ModifyCorrectionFast(double h, Teuchos::RCP<const Vector> res,
                       Teuchos::RCP<const Vector> u, Teuchos::RCP<Vector> du) override
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


};



#endif
