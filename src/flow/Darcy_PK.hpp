#ifndef __Darcy_PK_hpp__
#define __Darcy_PK_hpp__

#include "Teuchos_RCP.hpp"

#include "Epetra_Vector.h"
#include "AztecOO.h"

#include "Flow_PK.hpp"
#include "Flow_State.hpp"
#include "DarcyProblem.hpp"

namespace Amanzi
{

class Darcy_PK : public Flow_PK
{

public:
  Darcy_PK(Teuchos::ParameterList&, const Teuchos::RCP<const Flow_State>);

  ~Darcy_PK ();

  int advance_to_steady_state();
  void commit_state(Teuchos::RCP<Flow_State>) {}

  // After a successful advance() the following routines may be called.

  // Returns a reference to the cell pressure vector.
  const Epetra_Vector& Pressure() const { return *pressure; }

  // Returns a reference to the Darcy face flux vector.
  const Epetra_Vector& Flux() const { return *darcy_flux; }

  // Computes the components of the Darcy velocity on cells.
  void GetVelocity (Epetra_MultiVector &q) const
      { problem->DeriveDarcyVelocity(*solution, q); }

private:

  Teuchos::RCP<const Flow_State> FS;
  Teuchos::RCP<FlowBC> bc;

  DarcyProblem *problem;
  AztecOO *solver;

  Epetra_Vector *solution;   // full cell/face solution
  Epetra_Vector *pressure;   // cell pressures
  Epetra_Vector *darcy_flux; // Darcy face fluxes

  int max_itr;      // max number of linear solver iterations
  double err_tol;   // linear solver convergence error tolerance
};

}

#endif
