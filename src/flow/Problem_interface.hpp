#ifndef __Problem_interface__
#define __Problem_interface__


#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"

#include "Teuchos_ParameterList.hpp"

#include "Darcy_problem.hpp"



// this is the interface to the pde problem as NOX expects it

class Problem_interface      : public NOX::Epetra::Interface::Required,
                               public NOX::Epetra::Interface::Jacobian,
                               public NOX::Epetra::Interface::Preconditioner

{

public:
 
  //! Constructor
  Problem_interface( Darcy_problem * Problem) :
    Problem_(Problem) {};

  //! Destructor
  ~Problem_interface() 
  {
  }

  bool computeF(const Epetra_Vector & x, Epetra_Vector & f,
                NOX::Epetra::Interface::Required::FillType F )
  {
    Problem_->ComputeF(x,f);

    return true;
  }
  
  bool computeJacobian(const Epetra_Vector & x, Epetra_Operator & Jac)
  {
    // this should never be called when using for example JFNK
    return true;
  }

  bool computePrecMatrix(const Epetra_Vector & x, Epetra_RowMatrix & M) 
  {
    // compute the preconditioning matrix here
    
    return true;
  }  
  
  bool computePreconditioner(const Epetra_Vector & x, Epetra_Operator & O,
                 Teuchos::ParameterList* pl)
  {
    return true;
  }  

private:
  
  Darcy_problem * Problem_;
  
}; /* class Problem_interface */



#endif
