#ifndef _BDF2_DAE_HPP_
#define _BDF2_DAE_HPP_

// This class is based on Neil Carlson's BDF2_DAE module 
// that is part of LANL's Truchas code. 

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"

#include "Epetra_Vector.h"
#include "Epetra_BlockMap.h"

#include "NKA.H"

#include "BDF2_State.hpp"
#include "BDF2_fnBase.hpp"


namespace BDF2 {

  class Dae : public Teuchos::VerboseObject<Dae>,
	      public Teuchos::ParameterListAcceptor
  {
    
  public:
    
    // Create the BDF2 Dae solver object, the nonlinear problem must
    // be defined in a class that derives from the virtual base class
    // fnBase.								
    // The map is passed in, so that the BDF2 Dae stepper knows what
    // kind of Epetra_Vector is needs to work with.
    // The parameter list plist is checked for validity in the constructor.
    Dae(fnBase& fn_, Epetra_BlockMap& map_);

    // initializes the state of the BDF2 stepper
    void set_initial_state(const double t, const Epetra_Vector& x, const Epetra_Vector& xdot);

    // after a successful BDF2 step, this method is used to commit
    // the new solution to the solution history
    void commit_solution(const double h, const Epetra_Vector& u);

    // based on the past time step sizes and the current error, this 
    // method computes a new time step length
    void select_step_size(const std::vector<double>& dt, const double perr, double& h); 

    // in the first time step, there is no solution history, so we need to boot-strap
    void trap_step_one(double h, Epetra_Vector& u, double& hnext, int& errc);

    // computes a BDF2 step
    void bdf2_step_gen(double h, Epetra_Vector& u, double& hnext, int& errc, bool ectrl);

    // a simple interface to trap_step_one and bdf2_step_gen
    void bdf2_step_simple(double h, Epetra_Vector& u, double& hnext, int& errc, bool ectrl);

    // attempts to take a BDF2 step, interfaces to bdf2_step_simple, will 
    // try several times with the suggested reduced time step to take the step
    void bdf2_step(double& h, double hmin, int mtries, Epetra_Vector& u, double& hnext); 

    // the nonlinear solver (uses NKA)
    void solve_bce(double t, double h, Epetra_Vector& u0, Epetra_Vector& u, int& errc);

    // returns the most recent time
    const double most_recent_time() 
    {
      return state.uhist->most_recent_time();
    }

    // write statistics about the time step
    void write_bdf2_stepping_statistics();

    // Overridden from ParameterListAccpetor
    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const&) ;
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList() ;
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;


  private:

    double rmin; 
    double rmax; 
    double margin;

    State state;
    nka* fpa;
    fnBase& fn;

    Epetra_BlockMap& map;

    Teuchos::RCP<Teuchos::ParameterList> paramList_;

    // constants
    const static double RMIN = 0.25;
    const static double RMAX = 4.0;
    const static double MARGIN = 3.0;

  };

}

#endif // _BDF2_DAE_HPP_
