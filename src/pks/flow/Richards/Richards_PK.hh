#ifndef __RICHARDS_PK_HH__
#define __RICHARDS_PK_HH__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_Vector.h"

#include "Vector.hh"
#include "PK.hh"
#include "State.hh"
#include "RichardsProblem.hh"
#include "RichardsModelEvaluator.hh"
#include "BDF2_Dae.hpp"

namespace Amanzi {
class Richards_PK : public PK {

public:
  Richards_PK(Teuchos::ParameterList &plist, Teuchos::RCP<State> &S);

  void initialize(Teuchos::RCP<State> &S);

  double get_dt() { return h; }

  bool advance(double dt, const Teuchos::RCP<State> &S0,
               Teuchos::RCP<State> &S1, Teuchos::RCP<Vector> &solution);

  void commit_state(double dt, Teuchos::RCP<State> &S) {}

  // // again, are these necessary?
  // // Returns a reference to the cell pressure vector.
  // const Epetra_Vector& Pressure() const { return *pressure_cells; }

  // // Returns a reference to the Richards face flux vector.
  // const Epetra_Vector& Flux() const { return *richards_flux; }

  // // make these private?
  // // Computes the components of the Richards velocity on cells.
  // void GetVelocity(Epetra_MultiVector &q) const
  //     { problem->DeriveDarcyVelocity(*solution, q); }

  // // Computes the fluid saturation on cells.
  // void GetSaturation(Epetra_Vector &s) const;

private:
  void advance_to_steady_state();

  Teuchos::RCP<RichardsProblem> problem;
  Teuchos::RCP<RichardsModelEvaluator> RME;

  Teuchos::RCP<BDF2::Dae> time_stepper;

  Teuchos::RCP<Epetra_Vector> solution;   // full cell/face solution
  Teuchos::RCP<Epetra_Vector> pressure_cells;   // cell pressures
  Teuchos::RCP<Epetra_Vector> pressure_faces;

  Teuchos::RCP<Epetra_Vector> richards_flux; // Darcy face fluxes

  int max_itr;      // max number of linear solver iterations
  double err_tol;   // linear solver convergence error tolerance
  int precon_freq;  // preconditioner update frequency

  bool steady_state;
  double ss_t0, ss_t1; // steady state start/stop times
  double h0; // initial timestep size
  double height0; // hydrostatic height at t0

  double h, hnext;

  Teuchos::ParameterList plist;
};

} // close namespace Amanzi

#endif
