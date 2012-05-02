#ifndef _BDF2TIMEINTEGRATOR_HH_
#define _BDF2TIMEINTEGRATOR_HH_

// This class is based on Neil Carlson's BDF2_DAE module
// that is part of LANL's Truchas code.

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"

#include "NonlinarKrylovAccelerator.hh"

#include "BDFTimeIntegrator.hh"
#include "BDF2State.hh"
#include "BDFFnBase.hh"


namespace Amanzi {

class BDF2TimeIntegrator : public BDFTimeIntegrator,
                           public Teuchos::VerboseObject<BDF2TimeIntegrator>
{

 public:

  // Create the BDF2 Dae solver object, the nonlinear problem must
  // be defined in a class that derives from the virtual base class
  // BDFFnBase.
  explicit BDF2TimeIntegrator(Amanzi::BDFFnBase* fn,
                              Teuchos::RCP<Teuchos::ParameterList> const& paramList,
                              Teuchos::RCP<Amanzi::TreeVector> initvector);

  ~BDF2TimeIntegrator();

  // initializes the state of the BDF2 stepper
  void set_initial_state(const double t, Teuchos::RCP<Amanzi::TreeVector> u, 
			 Teuchos::RCP<Amanzi::TreeVector> udot);

  // after a successful BDF2 step, this method is used to commit
  // the new solution to the solution history
  void commit_solution(const double h, Teuchos::RCP<Amanzi::TreeVector> u);

  // attempts to take a BDF2 step, interfaces to bdf2_step_simple, will
  // try several times with the suggested reduced time step to take the step
  // returns the suggested next time step
  double time_step(const double h, Teuchos::RCP<Amanzi::TreeVector> u);

  void reset();

  // returns the most recent time
  double time();

 private:
  // read initialization parameters from a parameter list
  void readParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);

  // based on the past time step sizes and the current error, this
  // method computes a new time step length
  double select_step_size(const std::vector<double>& dt, const double perr);

  // in the first time step, there is no solution history, so we need to boot-strap,
  // returns the next time step
  // returns the next time step
  double trap_step_one(double h, Teuchos::RCP<Amanzi::TreeVector> u);

  // computes a BDF2 step
  // returns the next time step
  double bdf2_step_gen(double h, Teuchos::RCP<Amanzi::TreeVector> u, bool ectrl);

  // a simple interface to trap_step_one and bdf2_step_gen
  // returns the next time step
  double bdf2_step_simple(double h, Teuchos::RCP<Amanzi::TreeVector> u, bool ectrl);

  // attempts to take a BDF2 step, interfaces to bdf2_step_simple, will
  // try several times with the suggested reduced time step to take the step
  // returns the next time step
  double bdf2_step(double h, double hmin, Teuchos::RCP<Amanzi::TreeVector> u);

  // the nonlinear solver (uses NKA)
  void solve_bce(double t, double h, Teuchos::RCP<Amanzi::TreeVector> u0, 
		 Teuchos::RCP<Amanzi::TreeVector> u);

  // write statistics about the time step
  void write_bdf2_stepping_statistics();

 private:

  double rmin_;
  double rmax_;
  double margin_;

  int mtries_;
  double hmin_;

  BDF2State state_;
  NonlinearKrylovAccelerator* fpa_;
  SolutionHistory* sh_;
  BDFFnBase& fn_;

  Teuchos::RCP<Amanzi::TreeVector> initvector_;

  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  // constants
  const static double RMIN_;
  const static double RMAX_;
  const static double MARGIN_;

};

}

#endif // _BDF2_HH_
