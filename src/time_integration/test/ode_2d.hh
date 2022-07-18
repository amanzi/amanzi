#ifndef ODE_2D_BASE_HH_
#define ODE_2D_BASE_HH_

#include "FnBaseDefs.hh"
#include <vector>
#include <cmath>

#include "Teuchos_RCP.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_DataAccess.h"

#include "PartitionFnBase.hh"


namespace Amanzi
{
  /**
   * @brief A simple 2d test problem
   * 
   * @tparam Interface interface to inherite from
   *                   must provide the given overrides to be valid:
   *                    FunctionalTimeDerivativeFull, FunctionalTimeDerivativeFast, FunctionalTimeDerivativeSlow  
   */
  template<class Interface>
  class ode_2d : public Interface
  {
    using Vector = Epetra_Vector;
    using Matrix = Epetra_CrsMatrix;
  protected :
    
    //Coefficents
    Teuchos::RCP<Matrix> A_partition_;
    Teuchos::RCP<Matrix> A_fast_;
    Teuchos::RCP<Matrix> A_slow_;
    Teuchos::RCP<Matrix> A_;
    int matrix_n = 2;
    int matrix_partition_n = 4;

    

    int row_index[2] = {0,1};

    //Mappings
    Teuchos::RCP<Epetra_Map> matrix_map;
    Teuchos::RCP<Epetra_Map> matrix_partition_map;

    /**
     * @brief Used for block wise scaling
     * 
     * Given a vector of M_blockamount * N_Blockamount scaling values
     * Scale the matrix block wise on 
     * 
     * Didn't find a wise block scaling code in Trilinos or possible transformation
     * 
     * @param scalings 
     * @param j_eval 
     */
    void block_scaling(int M_blockamount, int N_blockamount,int M_blocksize, int N_blocksize,std::vector<double> scalings, Teuchos::RCP<Matrix> j_eval){
      for (auto i = 0; i < M_blockamount; i++)
      {
        for (auto j = 0; j < N_blockamount; j++)
        {
          for (auto k = 0; k < M_blocksize; k++)
          {
            for (auto p = 0; p < N_blocksize; p++)
            {
              if ((*j_eval)[i * M_blockamount + k][j*N_blockamount + p] != 0)
              {
                (*j_eval)[i * M_blockamount + k][j*N_blockamount + p] *= scalings[i * M_blockamount + j];
              }
              
            }
            
          }
          
        }
        
      }
      
    }

  public:
      ode_2d(Teuchos::RCP< Matrix> A_fast, Teuchos::RCP< Matrix> A_slow, Epetra_MpiComm* comm);
      void FunctionalTimeDerivativeFull(const double t, const Teuchos::RCP<Vector> u, Teuchos::RCP<Vector> u_eval) override;
      void FunctionalTimeDerivativeSlow(const double t, const Teuchos::RCP<Vector> u, Teuchos::RCP<Vector> f) override;
      void FunctionalTimeDerivativeFast(const double t, const Teuchos::RCP<Vector> u, Teuchos::RCP<Vector> f) override;
      void exact_rhs(double t, Teuchos::RCP<Vector> u_eval);
  };

  
  template<class Interface>
  ode_2d<Interface>::ode_2d(Teuchos::RCP< Matrix> A_fast, Teuchos::RCP< Matrix> A_slow, Epetra_MpiComm* comm):
          A_fast_(A_fast), A_slow_(A_slow)
  {
    
    matrix_map = Teuchos::rcp( new Epetra_Map(2, 0, *comm));
    matrix_partition_map = Teuchos::rcp( new Epetra_Map(4, 0, *comm));
    
    A_partition_ =  Teuchos::rcp(new Matrix(Copy, *matrix_partition_map, 4, true));
    
    A_ =  Teuchos::rcp(new Matrix(Copy, *matrix_map, 2, true));

    //Initalize the Full A
    *A_ = *A_fast_;

    int numEntries = 0;
    double* vals = new double [2];
    int* indices = new int[2];

    for (int i = 0; i < A_->NumMyRows(); i++)
    {
      A_slow_->ExtractGlobalRowCopy(i, 2, numEntries, vals, indices);
      A_->SumIntoGlobalValues(i, numEntries, vals, indices);
    }

    //Initalize the Partition Matrix
    for (int i = 0; i < matrix_n; i++)
    {
      A_fast_->ExtractGlobalRowCopy(i, 2, numEntries, vals, indices);
      A_partition_->InsertGlobalValues(i, numEntries, vals, indices);
      A_partition_->InsertGlobalValues(i + matrix_n, numEntries, vals, indices );

      A_slow_->ExtractGlobalRowCopy(i, 2, numEntries, vals, indices);
      indices[0] += matrix_n;
      indices[1] += matrix_n;
      A_partition_->InsertGlobalValues(i, numEntries, vals, indices);
      A_partition_->InsertGlobalValues(i + matrix_n, numEntries, vals, indices );
    }
    
    A_fast_->FillComplete();
    A_slow_->FillComplete();
    A_->FillComplete();
    A_partition_->FillComplete();

  }

  
/**
 * @brief Evaluate the full RHS of the problem
 * 
 * @param t 
 * @param u 
 * @param u_eval results
 */
  template<class Interface>
  void ode_2d<Interface>::FunctionalTimeDerivativeFull(const double t, const Teuchos::RCP<Vector> u, Teuchos::RCP<Vector> u_eval) {
    A_->Multiply(false, *u, *u_eval);
  }

/**
 * @brief Slow addiatively split term
 * 
 * @param t 
 * @param u 
 * @param u_eval 
 */
  template<class Interface>
  void ode_2d<Interface>::FunctionalTimeDerivativeSlow(const double t, const Teuchos::RCP<Vector> u, Teuchos::RCP<Vector> u_eval){
    A_slow_->Multiply(false, *u, *u_eval);
  }

  template<class Interface>
  void ode_2d<Interface>::FunctionalTimeDerivativeFast(const double t, const Teuchos::RCP<Vector> u, Teuchos::RCP<Vector> u_eval){
    A_fast_->Multiply(false, *u, *u_eval);
  }




/**
 * @brief Exact solution to the ODE
 * 
 * y' = A y
 * 
 * where A \in R^{2 x 2} & y \in R^2
 * 
 * @param t 
 * @param u_eval 
 */
  template<class Interface>
  void ode_2d<Interface>::exact_rhs(double t, Teuchos::RCP<Vector> u_eval){
    if ((*A_)[0][1] == 0 && (*A_)[1][0] == 0)
    {
      (*u_eval)[0] = exp((*A_)[0][0] * t);
      (*u_eval)[1] = exp((*A_)[1][1] * t);

      return;
    }
    
    
    double temp_e = exp((1.0/2.0) * t *  ((*A_)[0][0] + (*A_)[1][1]));
    
    double temp_sqrt = sqrt(4.0 * (*A_)[0][1] * (*A_)[1][0] + pow(((*A_)[0][0] - (*A_)[1][1]),2.0));
    
    double temp_cosh = temp_e * cosh((1.0/2.0) * t * temp_sqrt);

    double temp_sinh = temp_e * sinh((1.0/2.0) * t * temp_sqrt) / temp_sqrt;

    (*u_eval)[0] = temp_cosh + temp_sinh * ((*A_)[0][0] + 2.0*(*A_)[0][1] - (*A_)[1][1]);
    (*u_eval)[1] = temp_cosh + temp_sinh * (-(*A_)[0][0] + 2.0*(*A_)[1][0] + (*A_)[1][1]);
  }

}





#endif