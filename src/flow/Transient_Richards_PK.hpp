#ifndef __Transient_Richards_PK_hpp__
#define __Transient_Richards_PK_hpp__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_Vector.h"

#include "Flow_PK.hpp"
#include "Flow_State.hpp"
#include "RichardsProblem.hpp"

class Transient_Richards_PK : public Flow_PK
{

public:
  Transient_Richards_PK(Teuchos::ParameterList&, const Teuchos::RCP<const Flow_State>);

  ~Transient_Richards_PK ();

  int advance();
  void commit_state(Teuchos::RCP<Flow_State>) {}

  // After a successful advance() the following routines may be called.

  // Returns a reference to the cell pressure vector.
  const Epetra_Vector& Pressure() const { return *pressure; }

  // Returns a reference to the Richards face flux vector.
  const Epetra_Vector& RichardsFlux() const { return *richards_flux; }

  // Computes the components of the Richards velocity on cells.
  void GetRichardsVelocity(Epetra_MultiVector &q) const
      { problem->DeriveDarcyVelocity(*solution, q); }

  // Computes the fluid saturation on cells.
  void GetSaturation(Epetra_Vector &s) const;

private:

  Teuchos::RCP<const Flow_State> FS;
  Teuchos::RCP<FlowBC> bc;

  RichardsProblem *problem;

  Teuchos::RCP<Teuchos::ParameterList> nox_param_p;
  Teuchos::RCP<Teuchos::ParameterList> linsol_param_p;

  Epetra_Vector *solution;   // full cell/face solution
  Epetra_Vector *pressure;   // cell pressures
  Epetra_Vector *richards_flux; // Darcy face fluxes

  int max_itr;      // max number of linear solver iterations
  double err_tol;   // linear solver convergence error tolerance
  int precon_freq;  // preconditioner update frequency


};

#endif
