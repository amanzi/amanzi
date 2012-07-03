
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "BDF1TimeIntegrator.hh"
#include "BDFFnBase.hh"
#include "SolutionHistory.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"


namespace Amanzi {

BDF1TimeIntegrator::BDF1TimeIntegrator(BDFFnBase* fn,
                                       Teuchos::RCP<Teuchos::ParameterList> const& paramList,
                                       Teuchos::RCP<Amanzi::TreeVector> initvector) : fn_(*fn)  {
  // set the line prefix for output
  this->setLinePrefix("BDF1TimeIntegrator         ");

  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);

}

BDF1TimeIntegrator::~BDF1TimeIntegrator() {
  delete fpa_;
  delete sh_;
}


void BDF1TimeIntegrator::readParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList) {

  paramList_ = paramList;

  // read the parameter list and initialize the class
  state_.mitr = paramList_->get<int>("steady limit iterations");
  state_.maxitr = paramList_->get<int>("steady max iterations");
  state_.minitr = paramList_->get<int>("steady min iterations");
  state_.ntol = paramList_->get<double>("steady nonlinear tolerance");
  state_.hred = paramList_->get<double>("steady time step reduction factor");
  state_.hinc = paramList_->get<double>("steady time step increase factor");
  state_.hlimit = paramList_->get<double>("steady max time step");

  state_.maxpclag = paramList_->get<int>("steady max preconditioner lag iterations");

  // sanity check
  if ( ! ((state_.minitr < state_.maxitr) && (state_.maxitr < state_.mitr) )) {
    Errors::Message m("steady state paramters are wrong, we need min iterations < max iterations < limit iterations");
    Exceptions::amanzi_throw(m);
  }


  int maxv = state_.mitr-1;
  int mvec = 10;
  maxv = std::min<int>(maxv,mvec);

  double vtol = 0.05;
  fpa_ = new NonlinearKrylovAccelerator(maxv, vtol, *initvector_);

  // create the solution history object
  sh_ = new SolutionHistory(2, initvector_);
  state_.init_solution_history(sh_);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&*paramList_,this);

}


void BDF1TimeIntegrator::commit_solution(const double h, Teuchos::RCP<Amanzi::TreeVector> u) {

  double t = h + state_.uhist->most_recent_time();

  state_.uhist->record_solution(t, u);

  state_.hlast = h;
  state_.seq++;
  state_.freeze_count = std::max<int>(0, state_.freeze_count-1);

  state_.hmin = std::min<double>(h, state_.hmin);
  state_.hmax = std::max<double>(h, state_.hmax);

}


void BDF1TimeIntegrator::set_initial_state(const double t, Teuchos::RCP<Amanzi::TreeVector> x,
                                           Teuchos::RCP<Amanzi::TreeVector> xdot)  {
  state_.uhist->flush_history(t, x, xdot);
  state_.seq = 0;
  state_.usable_pc = false;
  state_.pclagcount = 0;
}



double BDF1TimeIntegrator::time_step(double h, Teuchos::RCP<Amanzi::TreeVector> u) {

  double tlast = state_.uhist->most_recent_time();
  double tnew = tlast + h;

  bool fresh_pc = false;

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  // start with guessing that the next time step should be the same as this one
  double hnext(h);

  // save the initial value so that we can restore the initial state
  // in case we need to bail
  Teuchos::RCP<Amanzi::TreeVector> usav = Teuchos::rcp(new Amanzi::TreeVector(std::string("u"),*u));
  *usav = *u;


  if (h < state_.hmin) {
    std::cout << h << " " << state_.hmin << std::endl;

    std::string msg = "BDF1 failed: Time step crash";
    Errors::Message m(msg);
    Exceptions::amanzi_throw(m);
  }

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
    *out << "BDF1 step " << state_.seq+1 << " T = "<< tlast << " H = " << h << std::endl;
  }

  // Predicted solution (initial value for the nonlinear solver)
  Teuchos::RCP<Amanzi::TreeVector> up = Teuchos::rcp(new Amanzi::TreeVector(std::string("up"),*u));

  if ( state_.uhist->history_size() > 1) {
    state_.uhist->interpolate_solution(tnew,  up);
  } else  {
    *up = *u;
  }

  // u at the start of the time step
  Teuchos::RCP<Amanzi::TreeVector> u0 = Teuchos::rcp(new Amanzi::TreeVector(std::string("u0"),*u));
  *u0 = *u;

  if (state_.pclagcount > state_.maxpclag) {
    state_.usable_pc = false;
    state_.pclagcount = 0;
  }

  if (!state_.usable_pc) {
    // update the preconditioner (we use the predicted solution at the end of the time step)
    state_.updpc_calls++;
    state_.pclagcount++;
    int errc = 0;
    try {
      fn_.update_precon (tnew, up, h);
    }
    catch (...) {
      std::string msg = "BDF1 failed: error while updating the preconditioner.";
      Errors::Message m(msg);
      Exceptions::amanzi_throw(m);
    }
    state_.usable_pc = true;
  }

  //  Solve the nonlinear BCE system.
  u = up; // Initial solution guess is the predictor.
  try {
    solve_bce ( tnew, h, u0, u);
  }
  catch (int itr) {
    // we end up in here either if the solver took too many iterations,
    // or if it took too few
    if (itr > state_.maxitr && itr <= state_.mitr) {
      hnext = std::min(state_.hred * h, state_.hlimit);
      state_.usable_pc = false;
    } else if (itr < state_.minitr) {
      hnext = std::min( state_.hinc * h, state_.hlimit);
      state_.usable_pc = false;
    } else if (itr > state_.mitr) {
      *u = *usav; // restore the original u
      state_.usable_pc = false;
      hnext = std::min(h, state_.hlimit);
      throw itr;
    }
  }

  return hnext;
}


void BDF1TimeIntegrator::solve_bce(double t, double h,  Teuchos::RCP<Amanzi::TreeVector> u0,
                                   Teuchos::RCP<Amanzi::TreeVector> u) {

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  fpa_->nka_restart();

  int itr = 0;

  Teuchos::RCP<Amanzi::TreeVector> du = Teuchos::rcp( new Amanzi::TreeVector(std::string("du"), *initvector_ ) );
  Teuchos::RCP<Amanzi::TreeVector> u_tmp = Teuchos::rcp( new Amanzi::TreeVector(std::string("u_tmp"), *initvector_ ) );

  Teuchos::RCP<Amanzi::TreeVector> preconditioned_f =
      Teuchos::rcp( new Amanzi::TreeVector(std::string("preconditioned_f"), *initvector_ ) );

  double error(0.0);
  double initial_error(0.0);

  do {

    // Check for too many nonlinear iterations.
    if (itr > state_.mitr) {
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "AIN BCE solve failed " << itr << " iterations (max), error = " << error << std::endl;
      }
      throw itr;
    }

    itr++;

    // Evaluate the preconditioned nonlinear function.
    state_.pcfun_calls++;

    // compute u_tmp = (u-u0)/h
    //*u_tmp = *u;
    //u_tmp->Update(-1.0/h,*u0,1.0/h);

    fn_.fun(t-h, t, u0, u, du);

    fn_.precon(du, u_tmp);

    // Accelerated correction.
    *preconditioned_f  = *u_tmp;        // copy preconditioned functional into appropriate data type

    fpa_->nka_correction(*du, *preconditioned_f);

    // Next solution iterate and error estimate.
    // FORTRAN:  u  = u - du
    u->Update(-1.0,*du,1.0);

    // Check the solution iterate for admissibility.
    if ( ! fn_.is_admissible(u) ) { // iterate is bad; bail.
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "AIN BCE solve FAILED: inadmissible solution iterate: itr=" << itr << std::endl;
      }
      throw std::string("solution iterate is inadmissible");
    }

    error = fn_.enorm(u, du);
    if (itr == 1) initial_error = error;

    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
      *out << itr << ": error = " << error << std::endl;
    }

    if (itr > 1 && error > initial_error*state_.elimit ) {
      // the solver threatening to diverge
      throw state_.mitr+1;
    }


    // Check for convergence
    if (error < state_.ntol)   {
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "AIN BCE solve succeeded: " << itr << " iterations, error = "<< error << std::endl;
      }

      if ((itr < state_.minitr) || (itr > state_.maxitr)) {
        throw itr;
      }
      return;
    }
  }
  while (true);

}


void BDF1TimeIntegrator::write_bdf1_stepping_statistics() {
  std::ostringstream oss;

  oss.flush();
  oss.setf(std::ios::scientific, std::ios::floatfield);

  oss << "STEP=";
  oss.fill('0');
  oss.width(5);
  oss << std::right << state_.seq;
  oss << " T=";
  oss.precision(5);
  oss.width(11);
  oss << std::left << state_.uhist->most_recent_time();
  oss << " H=";
  oss.precision(5);
  oss.width(11);
  oss << std::left << state_.hlast;

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true))  {
    *out << oss.str() << std::endl;
  }
}

}
