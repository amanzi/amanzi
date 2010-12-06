#ifndef __RICHARDSNOXINTERFACE__
#define __RICHARDSNOXINTERFACE__

#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"

#include "RichardsProblem.hpp"

// This is the interface to the problem as NOX expects it.

class RichardsNoxInterface : public NOX::Epetra::Interface::Required,
                             public NOX::Epetra::Interface::Jacobian,
                             public NOX::Epetra::Interface::Preconditioner
{

public:
 
  RichardsNoxInterface(RichardsProblem *problem) : problem_(problem) {}

  ~RichardsNoxInterface() {}

  // Compute F(x).  Virtual function from "Required".
  bool computeF(const Epetra_Vector &x, Epetra_Vector &f, FillType flag);
  
  // Compute the Jacobian DF(x).  Virtual function from "Jacobian"
  bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &J);

  // Compute the preconditioning operator.  Virtual function from "Preconditioner".
  bool computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M, Teuchos::ParameterList *params);

private:
  
  RichardsProblem *problem_;
  
};

#endif
