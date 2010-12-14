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
 
  RichardsNoxInterface(RichardsProblem *problem) : problem_(problem) { lag_prec = 1; lag_count=0; }

  ~RichardsNoxInterface() {}

  // Compute F(x).  Virtual function from "Required".
  bool computeF(const Epetra_Vector &x, Epetra_Vector &f, FillType flag);
  
  // Compute the Jacobian DF(x).  Virtual function from "Jacobian"
  bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &J);

  // Compute the preconditioning operator.  Virtual function from "Preconditioner".
  bool computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M, Teuchos::ParameterList *params);

  
  inline void setPrecLag(int lag_prec_) { lag_prec = lag_prec_;}
  inline void resetPrecLagCounter () { lag_count = 0; }
  inline int getPrecLag () { return lag_prec; }
  inline int getPrecLagCounter () { return lag_count; }

private:
  
  RichardsProblem *problem_;

  int lag_prec;  // the preconditioner is lagged this many times before it is recomputed
  int lag_count; // this counts how many times the preconditioner has been lagged
  
};

#endif
