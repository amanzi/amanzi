
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "ImplicitTIBDF2.hh"
#include "SolutionHistory.hh"
#include "ImplicitTIBDF2fnBase.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

namespace Amanzi {

const double ImplicitTIBDF2::RMIN_ = 0.25;
const double ImplicitTIBDF2::RMAX_ = 4.0;
const double ImplicitTIBDF2::MARGIN_ = 3.0;

ImplicitTIBDF2::ImplicitTIBDF2(ImplicitTIBDF2fnBase& fn, Teuchos::RCP<Amanzi::TreeVector> initvector) :
    fn_(fn) {
  initvector_ = Teuchos::rcp( new Amanzi::TreeVector(std::string("bdf2 initvector"), *initvector ) );
  
  // set the line prefix for output
  this->setLinePrefix("ImplicitTIBDF2         ");
  
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);
}

ImplicitTIBDF2::~ImplicitTIBDF2() {
  delete fpa_;
  delete sh_;
}

void ImplicitTIBDF2::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList) {
  TEST_FOR_EXCEPT(is_null(paramList));
  
  // Validate and set the parameter defaults.  Here, the parameters are
  // validated and the state of *this is not changed unless the parameter
  // validation succeeds.  Also, any validators that are defined for various
  // parameters are passed along so that they can be used in extracting
  // values!

  paramList->validateParametersAndSetDefaults(*this->getValidParameters(),0);
  paramList_ = paramList;

  // make sure that the parameter list is actually valid (this is probably redundant)
  paramList_->validateParameters(*this->getValidParameters());

  // read the parameter list and initialize the class
  state_.mitr = paramList_->get<int>("Nonlinear solver max iterations");
  state_.ntol = paramList_->get<double>("Nonlinear solver tolerance");

  mtries_ = paramList_->get<int>("Maximum number of BDF tries",20);

  int maxv = state_.mitr-1;
  int mvec = paramList_->get<int>("NKA max vectors");
  maxv = std::min<int>(maxv,mvec);

  // Initialize the FPA structure.
  double vtol = paramList_->get<double>("NKA drop tolerance");
  fpa_ = new NonlinearKrylovAccelerator(maxv, vtol, *initvector_);

  // create the solution history object
  sh_ = new SolutionHistory(3, initvector_);
  state_.init_solution_history(sh_);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
}



Teuchos::RCP<Teuchos::ParameterList> ImplicitTIBDF2::getNonconstParameterList() {
  return paramList_;
}



Teuchos::RCP<Teuchos::ParameterList> ImplicitTIBDF2::unsetParameterList() {
  Teuchos::RCP<Teuchos::ParameterList> paramList = paramList_;
  paramList_ = Teuchos::null;
  return paramList;
}



Teuchos::RCP<const Teuchos::ParameterList> ImplicitTIBDF2::getParameterList() const {
  return paramList_;
}



Teuchos::RCP<const Teuchos::ParameterList> ImplicitTIBDF2::getValidParameters() const {
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;
  static RCP<const ParameterList> validParams;
  if (is_null(validParams)) {
    RCP<ParameterList>
        pl = Teuchos::rcp(new ParameterList("Time Stepping"));
    Teuchos::setIntParameter("Nonlinear solver max iterations",
                             10,
                             "The maximum number of nonlinear iterations per invocation of the nonlinear solver.",
                             &*pl);
    Teuchos::setIntParameter("NKA max vectors",
                             5,
                             "The size of the vectorspace for the nonlinear Krylov accelerator.",
                             &*pl);
    Teuchos::setDoubleParameter("Nonlinear solver tolerance",
                                0.001,
                                "The tolerance for the nonlinear solver.",
                                &*pl);
    Teuchos::setDoubleParameter("NKA drop tolerance",
                                0.05,
                                "Drop tolerance for the nonlinear Krylov accelerator.",
                                &*pl);
    Teuchos::setIntParameter("Maximum number of BDF tries",
                             20,
                             "The number of failed BDF attempts we are tolerating.",
                             &*pl);
    Teuchos::setupVerboseObjectSublist(&*pl);
    validParams = pl;
  }
  return validParams;
}



void ImplicitTIBDF2::commit_solution(const double h, Teuchos::RCP<Amanzi::TreeVector> u) {
  
  double t = h + state_.uhist->most_recent_time();
  
  state_.uhist->record_solution(t, u);

  state_.hlast = h;
  state_.seq++;
  state_.freeze_count = std::max<int>(0, state_.freeze_count-1);

  state_.hmin = std::min<double>(h, state_.hmin);
  state_.hmax = std::max<double>(h, state_.hmax);
}




void ImplicitTIBDF2::select_step_size(const std::vector<double>& dt, const double perr, double& h) {
  const double tol = 0.001;

  ASSERT(dt.size() == 3);

  double a = 0.5*dt[0]*dt[1]*dt[2]/std::max<double>(perr,0.001);
  h = dt[0];

  double dh;
  do {
    double phi  = h*(h + dt[0])*(h + dt[1]) - a;
    double dphi = (2.0*h + dt[0])*(h + dt[1]) + h*(h + dt[0]);
    
    dh = phi / dphi;
    h = h - dh;
  } while (abs(dh) / h >= tol);  
}




void ImplicitTIBDF2::set_initial_state(const double t, Teuchos::RCP<Amanzi::TreeVector> x,
                                       Teuchos::RCP<Amanzi::TreeVector> xdot)  {
  // ASSERT(x.Map().PointSameAs( xdot.Map() ) );
  
  state_.uhist->flush_history(t, x, xdot);
  state_.seq = 0;
}


void ImplicitTIBDF2::bdf2_step(double& h, double hmin, Teuchos::RCP<Amanzi::TreeVector> u, double& hnext) {
  
  ASSERT(hmin>=0.0);
  ASSERT(hmin<=h);
  ASSERT(mtries_>=1);

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  int tries = 0;
  do {
    tries++;
    
    // Check for too many attempts at a single step.
    if (tries > mtries_) {
      Errors::Message m("BDF2 step failed");
      Exceptions::amanzi_throw(m);
    }
    
    // Check for a too-small step size.
    if (h < hmin) {
      Errors::Message m("BDF2 step size too small");
      Exceptions::amanzi_throw(m);
    }
    
    // Attempt a BDF2 step.
    int errc = 0;
    bdf2_step_simple (h, u, hnext, errc, true);
    if (errc == 0) return;

    // Step failed; try again with the suggested step size.
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
      *out << "Changing H by a factor of " << hnext/h << std::endl;
    }
    h = hnext;
  } while (true);
}



void ImplicitTIBDF2::bdf2_step_simple(double h, Teuchos::RCP<Amanzi::TreeVector> u,
                                      double& hnext, int& errc, bool ectrl) {
  ASSERT(state_.seq>=0);
  ASSERT(h>0);

  double mitr_sav;

  switch (state_.seq) {
    case 0:
      trap_step_one(h, u, hnext, errc);
      break;
      
    case 1:
    case 2:
      bdf2_step_gen(h, u, hnext, errc, false);
      break;
      
    default:
      bdf2_step_gen(h, u, hnext, errc, ectrl);
      break;
  }
}


void ImplicitTIBDF2::trap_step_one(double h, Teuchos::RCP<Amanzi::TreeVector> u, double& hnext, int& errc) {
  
  double tlast = state_.uhist->most_recent_time();
  double t = tlast + h;
  double etah = 0.5 * h;
  double t0 = tlast + etah;

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
    *out << "Trap step " << state_.seq+1 << " tlast = " << tlast << " h = " << h << " etah = " << etah << std::endl;
  }
  
  
  // Predicted solution and base point for the BCE step.
  Teuchos::RCP<Amanzi::TreeVector> u0 = Teuchos::rcp( new Amanzi::TreeVector(std::string("u0"), *initvector_ ) );

  state_.uhist->interpolate_solution (t,  u, 1);
  state_.uhist->interpolate_solution (t0, u0, 1);

  // Update the preconditioner.
  state_.updpc_calls++;
  fn_.update_precon (t, u, etah, errc);
  if (errc != 0) {
    state_.updpc_failed++;
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
      *out << "Preconditioner update failed at T = " << t << "  ETAH = " << etah << std::endl;
    }
    hnext = 0.1 * h; // want to quickly find an acceptably small step size.
    errc = -1;
    return;
  }
  
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
    *out << "Preconditioner updated at T = " << t << " etah = " << etah << std::endl;
  }
  state_.hpc = etah;
  state_.usable_pc = true;
  
  // Solve the nonlinear BCE system.
  solve_bce ( t, etah, u0, u, errc);
  if (errc != 0) {
    state_.failed_bce++;
    hnext = 0.1 * h;  // want to quickly find an acceptably small step size.
    state_.freeze_count = 1;
    errc = 1;
    return;
  }
  
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
    *out << "Step accepted: no local error control" << std::endl;
  }
  state_.freeze_count = 2;
  hnext = h;
  errc = 0;
}


void ImplicitTIBDF2::bdf1_step(double h, Teuchos::RCP<Amanzi::TreeVector> u) {
  double tlast = state_.uhist->most_recent_time(); 
  double t = tlast + h;
  int errc;

  Teuchos::RCP<Amanzi::TreeVector> up = Teuchos::rcp( new Amanzi::TreeVector(std::string("up"), *initvector_ ) );
  Teuchos::RCP<Amanzi::TreeVector> u0 = Teuchos::rcp( new Amanzi::TreeVector(std::string("u0"), *initvector_ ) );

  state_.uhist->interpolate_solution(t,  up, 2);  
  state_.uhist->interpolate_solution(tlast, u0, 1);

  fn_.update_precon (t, up, h, errc);
  *u = *up; // Initial solution guess is the predictor.
  
  solve_bce ( t, h, u0, u, errc );
  
}




void ImplicitTIBDF2::bdf2_step_gen(double h, Teuchos::RCP<Amanzi::TreeVector> u, double& hnext, int& errc, bool ectrl) {
  
  double tlast = state_.uhist->most_recent_time();
  double t = tlast + h;
  double eta = (state_.hlast + h) / (state_.hlast + 2.0 * h);
  double etah = eta * h;
  double t0 = tlast + (1.0 - eta)*h;

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
    *out << "BDF2 step " << state_.seq+1 << " T = "<< tlast << " H = " << h << " ETAH = " << etah << std::endl;
  }
  bool fresh_pc = false;
  
  // If the PC step size is too different than the current step size we tag
  // it as unusable in order to preempt a possible nonlinear solve failure.
  if (state_.usable_pc) {
    if (state_.hpc/etah > 1.0 + MARGIN_) state_.usable_pc = false;
    if (etah/state_.hpc > 1.0 + MARGIN_) state_.usable_pc = false;
  }
  
  // Predicted solution and base point for BCE step.
  Teuchos::RCP<Amanzi::TreeVector> up = Teuchos::rcp( new Amanzi::TreeVector(std::string("up"), *initvector_ ) );
  Teuchos::RCP<Amanzi::TreeVector> u0 = Teuchos::rcp( new Amanzi::TreeVector(std::string("u0"), *initvector_ ) );

  state_.uhist->interpolate_solution(t,  up, 2);
  state_.uhist->interpolate_solution(t0, u0, 1);

  // check the predicted solution for admissibility
  // and update preconditioner/halve time step if not
  if ( ! fn_.is_admissible(up) ) {
    state_.updpc_calls++;
    fn_.update_precon(t, up, etah, errc);
    if ( errc != 0) { // update failed 
      state_.updpc_failed++;
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "Preconditioner update FAILED at T=" << t << ", ETAH=" << etah << std::endl;
      }
      hnext = 0.5 * h;
      
      errc = -1;
      return;
    }
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
      *out << "Preconditioner updated at T=" << t << ", ETAH=" << etah << std::endl;
    }
    
    state_.hpc = etah;
    state_.usable_pc = true;
    fresh_pc = true;
  }
  
  // BCE loop:
  do {
    // Update the preconditioner if necessary.
    if (! state_.usable_pc) {
      state_.updpc_calls++;
      fn_.update_precon (t, up, etah, errc);
      if (errc != 0) { // update failed; cut h and return error condition
        state_.updpc_failed++;
        if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
          *out << "Preconditioner update failed, T = "<< t << " etah = " << etah << std::endl;
        }
        hnext = 0.5 * h;
        errc = -1;
        return;
      }
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "Preconditioner updated at T = "<< t << " etah = " << etah <<std::endl;
      }
      state_.hpc = etah;
      state_.usable_pc = true;
      fresh_pc = true;
    }
    
    //  Solve the nonlinear BCE system.
    *u = *up; // Initial solution guess is the predictor.
    solve_bce ( t, etah, u0, u, errc );
    if (errc == 0) break;  // leave the do-while loop, the BCE step was successful.
    
    if (fresh_pc) { // preconditioner was fresh; cut h and return error condition.
      state_.failed_bce++;
      hnext = 0.5 * h;
      state_.freeze_count = 1;
      errc = 1;
      return;
    } else { // update the preconditioner and retry the nonlinear solve.
      state_.retried_bce++;
      state_.usable_pc = false;
      continue; // cycle
    }
  } while (true);
  

  //bool predictor_error = true;
  bool predictor_error = ectrl;
  
  if (predictor_error) {
    
    // Predictor error control.
    // FORTRAN: u0 = u - up
    *u0 = *u;
    u0->Update(-1.0, *up, 1.0);
    
    
    double perr = fn_.enorm(u, u0);
    if (perr < 4.0) { // accept the step.
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "Step accepted, perr = " << perr <<std::endl;
      }
      errc = 0;
    } else { // reject the step; cut h and return error condition.
      state_.rejected_steps++;
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "Step REJECTED, perr = " << perr << std::endl;
      }
      hnext = 0.5 * h;
      state_.freeze_count = 1;
      errc = 2;
      return;
    }
    
    // Select the next step size based on the predictor error and past step
    // size history, but don't change the step size by too great a factor.
    // FORTRAN:  dt(1) = h
    // FORTRAN:  dt(2:) = h + time_deltas(this%uhist)
    
    std::vector<double> dt(3);
    std::vector<double> dt0(2); // to store the time deltas from the solution history
    state_.uhist->time_deltas(dt0);
    
    dt[0] = h;
    dt[1] = h + dt0[0];
    dt[2] = h + dt0[1];

    select_step_size (dt, perr, hnext);

    hnext = std::max<double>(RMIN_*h, std::min<double>(RMAX_*h, hnext));
    if (state_.freeze_count != 0) hnext = std::min<double>(h, hnext);
  } else {
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
      *out << "Step accepted: no local error control" << std::endl;
    }
    hnext = h;
    errc = 0; 
  }
}


void ImplicitTIBDF2::solve_bce(double t, double h, Teuchos::RCP<Amanzi::TreeVector> u0,
                               Teuchos::RCP<Amanzi::TreeVector> u, int& errc) {
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  fpa_->nka_restart();

  int itr = 0;

  Teuchos::RCP<Amanzi::TreeVector> du = Teuchos::rcp( new Amanzi::TreeVector(std::string("du"), *initvector_ ) );
  Teuchos::RCP<Amanzi::TreeVector> u_tmp = Teuchos::rcp( new Amanzi::TreeVector(std::string("u_tmp"), *initvector_ ) );


  //  Teuchos::RCP<NOX::Epetra::Vector> preconditioned_f =
  //    Teuchos::rcp(new NOX::Epetra::Vector(du, NOX::ShapeCopy));

  do {
    double error;
    
    // Check for too many nonlinear iterations.
    if (itr >= state_.mitr)  {
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "AIN BCE solve failed " << itr << " iterations (max), error = " << error << std::endl;
      }
      errc = 1;
      return;
    }

    itr++;

    // Evaluate the preconditioned nonlinear function.
    state_.pcfun_calls++;

    // compute u_tmp = (u-u0)/h
    *u_tmp = *u;
    u_tmp->Update(-1.0/h,*u0,1.0/h);

    fn_.fun(t, u, u_tmp, du);

    fn_.precon(du, u_tmp);

    fpa_->nka_correction(*u_tmp, *du);
    
    // Next solution iterate and error estimate.
    // FORTRAN:  u  = u - du
    u->Update(-1.0,*du,1.0);

    // Check the solution iterate for admissibility.
    if ( ! fn_.is_admissible(u) ) { // iterate is bad; bail.
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "AIN BCE solve FAILED: inadmissible solution iterate: itr=" << itr << std::endl;
        errc = 2;
        return;
      }
    }
    
    error = fn_.enorm(u, du);
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
      *out << itr << ": error = " << error << std::endl;
    }
    
    // Check for convergence.
    if (((error < state_.ntol) && (itr > 1)) || (error < 0.01 * state_.ntol)) {
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "AIN BCE solve succeeded: " << itr << " iterations, error = "<< error << std::endl;
      }
      errc = 0;
      return;
    }
  } while (true);

}


void ImplicitTIBDF2::write_bdf2_stepping_statistics() {

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
  oss << " NFUN:NPC=";
  oss.fill('0');
  oss.width(4);
  oss << std::right << state_.pcfun_calls;
  oss << ":";
  oss.fill('0');
  oss.width(4);
  oss << std::right << state_.updpc_calls;
  oss << " NPCF:NNR:NNF:NSR=";
  oss.fill('0');
  oss.width(4);
  oss << std::right << state_.updpc_failed;
  oss << ":";
  oss.fill('0');
  oss.width(2);
  oss << std::right << state_.retried_bce;
  oss << ":";
  oss.fill('0');
  oss.width(2);
  oss << std::right << state_.failed_bce;
  oss << ":";
  oss.fill('0');
  oss.width(2);
  oss << std::right << state_.rejected_steps;


  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true)) {
    *out << oss.str() << std::endl;
  }
}

}
