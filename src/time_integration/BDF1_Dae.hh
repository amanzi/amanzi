#ifndef _BDF1_DAE_HH_
#define _BDF1_DAE_HH_

// This class is based on Neil Carlson's BDF2_DAE module
// that is part of LANL's Truchas code.

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"

#include "Epetra_Vector.h"
#include "Epetra_BlockMap.h"

#include "NKA.H"
#include "BDF1_State.hh"
#include "BDF2_fnBase.hpp"


namespace Amanzi {

class BDF1Dae : public Teuchos::VerboseObject<BDF1Dae>,
                public Teuchos::ParameterListAcceptor {
 public:
  // Create the BDF1Dae solver object, the nonlinear problem must
  // be defined in a class that derives from the virtual base class
  // fnBase.
  // The map is passed in, so that the BDF2 Dae stepper knows what
  // kind of Epetra_Vector is needs to work with.
  // The parameter list plist is checked for validity in the constructor.
  BDF1Dae(BDF2::fnBase& fn_, const Epetra_BlockMap& map_);
  ~BDF1Dae();

  // initializes the state of the BDF2 stepper
  void set_initial_state(const double t, const Epetra_Vector& x, const Epetra_Vector& xdot);

  // after a successful BDF2 step, this method is used to commit
  // the new solution to the solution history
  void commit_solution(const double h, const Epetra_Vector& u);

  // based on the past time step sizes and the current error, this
  // method computes a new time step length
  void select_step_size(const std::vector<double>& dt, const double perr, double& h);

  // computes a BDF1 step
  void bdf1_step(double h, Epetra_Vector& u, double& hnext);

  // the nonlinear solver (uses NKA)
  void solve_bce(double t, double h, Epetra_Vector& u0, Epetra_Vector& u);

  // returns the most recent time
  const double most_recent_time() { return state.uhist->most_recent_time(); }

  // write statistics about the time step
  void write_bdf1_stepping_statistics();

  // Overridden from ParameterListAccpetor
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const&);
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

 private:
  double rmin;
  double rmax;
  double margin;

  int mtries;

  BDF1State state;
  nka* fpa;
  BDF2::SolutionHistory* sh_;
  BDF2::fnBase& fn;

  const Epetra_BlockMap& map;

  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  // constants
  const static double RMIN;
  const static double RMAX;
  const static double MARGIN;
};

}

#endif // _BDF2_DAE_HPP_
