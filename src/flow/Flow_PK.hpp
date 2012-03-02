#ifndef __Flow_PK_hpp__
#define __Flow_PK_hpp__

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Flow_State.hpp"

namespace Amanzi
{

// Pure virtual base class that the Darcy and Richards flow
// process kernels should derive from.

class Flow_PK {
public:
  virtual int advance_to_steady_state() = 0;
  
  virtual int advance_transient(double h) = 0; 
  virtual int advance_steady(double h) = 0; 
  virtual int init_transient(double t0, double h0) = 0; 
  virtual int init_steady(double t0, double h0) = 0; 
  
  virtual void commit_state(Teuchos::RCP<Flow_State>, double ) = 0;
  virtual void commit_new_saturation(Teuchos::RCP<Flow_State>) = 0;

  // After a successful advance() the following routines may be called.

  // Returns a reference to the cell pressure vector.
  virtual const Epetra_Vector& Pressure() const = 0;

  // Returns a reference to the Darcy face flux vector.
  virtual const Epetra_Vector& Flux() const = 0;

  // Computes the components of the Darcy velocity on cells.
  virtual void GetVelocity(Epetra_MultiVector &q) const = 0;

  virtual double get_flow_dT() = 0; 
  
};

}

#endif
