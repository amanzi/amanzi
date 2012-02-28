#ifndef _BDF1TIMEINTEGRATOR_HH_
#define _BDF1TIMEINTEGRATOR_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"

#include "NonlinarKrylovAccelerator.hh"
#include "SolutionHistory.hh"

#include "BDFTimeIntegrator.hh"
#include "BDF1State.hh"
#include "BDFFnBase.hh"


namespace Amanzi {

class BDF1TimeIntegrator : public BDFTimeIntegrator,
                           public Teuchos::VerboseObject<BDF1TimeIntegrator>

{

 public:

  // Create the BDF2 Dae solver object, the nonlinear problem must
  // be defined in a class that derives from the virtual base class
  // FnBase.

  BDF1TimeIntegrator(BDFFnBase* fn, 
                     Teuchos::RCP<Teuchos::ParameterList> const& paramList,
                     Teuchos::RCP<Amanzi::TreeVector> initvector);

  ~BDF1TimeIntegrator();

  // initializes the state of the BDF2 stepper
  void set_initial_state(const double h, Teuchos::RCP<Amanzi::TreeVector> x, 
                         Teuchos::RCP<Amanzi::TreeVector> xdot);

  // after a successful BDF2 step, this method is used to commit
  // the new solution to the solution history
  void commit_solution(const double h, Teuchos::RCP<Amanzi::TreeVector> u);

  // computes a BDF1 step
  double time_step(double h, Teuchos::RCP<Amanzi::TreeVector> x);  

  // this method will reset the memory of the time
  // time integrator
  void reset();

  // returns the most recent time
  double time();
  
  
 private:

  // the nonlinear solver (uses NKA)
  void solve_bce(double t, double h,  Teuchos::RCP<Amanzi::TreeVector> u0, Teuchos::RCP<Amanzi::TreeVector> u);

  // initialize from a paramter list
  void readParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);

  // write statistics about the time step
  void write_bdf1_stepping_statistics();

  int mtries_;
  double hmin_;

  BDF1State state_;
  NonlinearKrylovAccelerator* fpa_;
  SolutionHistory* sh_;
  BDFFnBase& fn_;

  Teuchos::RCP<Amanzi::TreeVector> initvector_;

  Teuchos::RCP<Teuchos::ParameterList> paramList_;

};

}

#endif // _BDF1TIMEINTEGRATOR_HH_
