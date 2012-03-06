#ifndef _IMPLICITTIBDF2_HH_
#define _IMPLICITTIBDF2_HH_

// This class is based on Neil Carlson's BDF2_DAE module 
// that is part of LANL's Truchas code. 

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"

#include "Epetra_Vector.h"
#include "Epetra_BlockMap.h"

#include "NonlinarKrylovAccelerator.hh"

#include "ImplicitTIBDF2State.hh"
#include "ImplicitTIBDF2fnBase.hh"


namespace Amanzi {

  class ImplicitTIBDF2 : public Teuchos::VerboseObject<ImplicitTIBDF2>,
			 public Teuchos::ParameterListAcceptor
  {
    
  public:
    
    // Create the BDF2 Dae solver object, the nonlinear problem must
    // be defined in a class that derives from the virtual base class
    // fnBase.								
    // The map is passed in, so that the BDF2 Dae stepper knows what
    // kind of Epetra_Vector is needs to work with.
    // The parameter list plist is checked for validity in the constructor.
    ImplicitTIBDF2(ImplicitTIBDF2fnBase& , Teuchos::RCP<Amanzi::TreeVector> );
    
    ~ImplicitTIBDF2();

    // initializes the state of the BDF2 stepper
    void set_initial_state(const double t, Teuchos::RCP<Amanzi::TreeVector> x, Teuchos::RCP<Amanzi::TreeVector> xdot);

    // after a successful BDF2 step, this method is used to commit
    // the new solution to the solution history
    void commit_solution(const double h, Teuchos::RCP<Amanzi::TreeVector> u);

    // based on the past time step sizes and the current error, this 
    // method computes a new time step length
    void select_step_size(const std::vector<double>& dt, const double perr, double& h); 

    // in the first time step, there is no solution history, so we need to boot-strap
    void trap_step_one(double h, Teuchos::RCP<Amanzi::TreeVector> u, double& hnext, int& errc);

    // computes a BDF2 step
    void bdf2_step_gen(double h, Teuchos::RCP<Amanzi::TreeVector> u, double& hnext, int& errc, bool ectrl);

    // a simple interface to trap_step_one and bdf2_step_gen
    void bdf2_step_simple(double h, Teuchos::RCP<Amanzi::TreeVector> u, double& hnext, int& errc, bool ectrl);

    // attempts to take a BDF2 step, interfaces to bdf2_step_simple, will 
    // try several times with the suggested reduced time step to take the step
    void bdf2_step(double& h, double hmin, Teuchos::RCP<Amanzi::TreeVector> u, double& hnext); 

    // takes a BDF1 (backward Euler) step of step size h)
    void bdf1_step(double h, Teuchos::RCP<Amanzi::TreeVector> u);
    // the nonlinear solver (uses NKA)
    void solve_bce(double t, double h, Teuchos::RCP<Amanzi::TreeVector> u0, Teuchos::RCP<Amanzi::TreeVector> u, int& errc);

    // returns the most recent time
    const double most_recent_time() 
    {
      return state_.uhist->most_recent_time();
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

    double rmin_; 
    double rmax_; 
    double margin_;
    
    int mtries_;

    ImplicitTIBDF2State state_;
    NonlinearKrylovAccelerator* fpa_;
    SolutionHistory* sh_;
    ImplicitTIBDF2fnBase& fn_;

    Teuchos::RCP<Amanzi::TreeVector> initvector_;

    Teuchos::RCP<Teuchos::ParameterList> paramList_;

    // constants
    const static double RMIN_;
    const static double RMAX_;
    const static double MARGIN_;

  };

}

#endif // _IMPLICITTIBDF2_HH_
