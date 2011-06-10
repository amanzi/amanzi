#ifndef __Transient_Richards_PK_hpp__
#define __Transient_Richards_PK_hpp__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_Vector.h"

#include "Flow_PK.hpp"
#include "Flow_State.hpp"
#include "RichardsProblem.hpp"
#include "RichardsModelEvaluator.hpp"
#include "BDF2_Dae.hpp"


class Transient_Richards_PK : public Flow_PK
{

public:
  Transient_Richards_PK(Teuchos::ParameterList&, const Teuchos::RCP<const Flow_State>);

  ~Transient_Richards_PK ();

  int advance_to_steady_state();
  int advance_transient();
  void commit_state(Teuchos::RCP<Flow_State>) {}

  // After a successful advance() the following routines may be called.

  // Returns a reference to the cell pressure vector.
  const Epetra_Vector& Pressure() const { return *pressure_cells; }

  // Returns a reference to the Richards face flux vector.
  const Epetra_Vector& Flux() const { return *richards_flux; }

  // Computes the components of the Richards velocity on cells.
  void GetVelocity(Epetra_MultiVector &q) const
      { problem->DeriveDarcyVelocity(*solution, q); }

  // Computes the fluid saturation on cells.
  void GetSaturation(Epetra_Vector &s) const;
  
  double get_flow_dT() { return h; }


private:

  Teuchos::RCP<const Flow_State> FS;
  Teuchos::RCP<FlowBC> bc;
  Teuchos::ParameterList &richards_plist;
  
  RichardsProblem *problem;
  RichardsModelEvaluator *RME;
  
  BDF2::Dae *time_stepper;

  Epetra_Vector *solution;   // full cell/face solution
  Epetra_Vector *pressure_cells;   // cell pressures
  Epetra_Vector *pressure_faces;
  Epetra_Vector *richards_flux; // Darcy face fluxes

  int max_itr;      // max number of linear solver iterations
  double err_tol;   // linear solver convergence error tolerance
  int precon_freq;  // preconditioner update frequency

  double ss_t0, ss_t1, ss_h0, ss_z;

  double h, hnext;

};

#endif
