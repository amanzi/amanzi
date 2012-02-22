#ifndef __Transient_Richards_PK_hpp__
#define __Transient_Richards_PK_hpp__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"

#include "Epetra_Vector.h"

#include "Flow_PK.hpp"
#include "Flow_State.hpp"
#include "RichardsProblem.hpp"
#include "RichardsModelEvaluator.hpp"
#include "BDF2_Dae.hpp"
#include "BDF1_Dae.hh"

namespace Amanzi {

class Transient_Richards_PK : public Flow_PK, public Teuchos::VerboseObject<Transient_Richards_PK> {
 public:
  Transient_Richards_PK(Teuchos::ParameterList&, const Teuchos::RCP<Flow_State>);
  ~Transient_Richards_PK ();

  int advance_to_steady_state();
  int advance_transient(double h);
  int advance_steady(double h);
  int init_transient(double t0, double h0);
  int init_steady(double t0, double h0);
  void commit_state(Teuchos::RCP<Flow_State>); 

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
  
  double get_flow_dT() { return hnext; }

private:
  void approximate_face_pressure(const Epetra_Vector& cell_pressure, Epetra_Vector& face_pressure);

  Teuchos::RCP<Flow_State> FS;
  Teuchos::ParameterList &plist;
  
  RichardsProblem *problem;
  RichardsModelEvaluator *RME;
  
  BDF2::Dae* time_stepper;
  BDF1Dae* steady_time_stepper;
  
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

}  // namespace Amanzi

#endif
