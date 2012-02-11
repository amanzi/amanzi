
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "NOX_Epetra_Vector.H"

#include "BDF1_Dae.hh"
#include "BDF2_SolutionHistory.hpp"
#include "BDF2_fnBase.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"


namespace Amanzi {

BDF1Dae::BDF1Dae(BDF2::fnBase& fn_, const Epetra_BlockMap& map_) :
    fn(fn_), map(map_)  {
  // set the line prefix for output
  this->setLinePrefix("BDF1Dae             ");

  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);
  
}

BDF1Dae::~BDF1Dae() {
  delete fpa;
  delete sh_;
}


void BDF1Dae::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList) {

  TEST_FOR_EXCEPT(is_null(paramList));

  // Validate and set the parameter defaults.  Here, the parameters are
  // validated and the state of *this is not changed unless the parameter
  // validation succeeds.  Also, any validators that are defined for various
  // parameters are passed along so theu they can be used in extracting
  // values!

  paramList->validateParametersAndSetDefaults(*this->getValidParameters(),0);
  paramList_ = paramList;

  // make sure that the parameter list is actually valid (this is probably redundant)
  paramList_->validateParameters(*this->getValidParameters());

  // read the parameter list and initialize the class
  state.mitr = paramList_->get<int>("steady limit iterations");
  state.maxitr = paramList_->get<int>("steady max iterations"); 
  state.minitr = paramList_->get<int>("steady min iterations");
  state.ntol = paramList_->get<double>("steady nonlinear tolerance");
  state.hred = paramList_->get<double>("steady time step reduction factor");
  state.hinc = paramList_->get<double>("steady time step increase factor");
  state.hlimit = paramList_->get<double>("steady max time step");
  state.maxpclag = paramList_->get<int>("steady max preconditioner lag iterations");

  // sanity check
  if ( ! ((state.minitr < state.maxitr) && (state.maxitr < state.mitr) )) {
    Errors::Message m("steady state paramters are wrong, we need min iterations < max iterations < limit iterations");
    Exceptions::amanzi_throw(m);        
  }

  
  int maxv = state.mitr-1;
  int mvec = 10;
  maxv = std::min<int>(maxv,mvec);

  // Initialize the FPA structure.
  // first create a NOX::Epetra::Vector to initialize nka with
  NOX::Epetra::Vector init_vector( Epetra_Vector(map), NOX::ShapeCopy );
  double vtol = 0.05;
  fpa = new nka(maxv, vtol, init_vector);

  // create the solution history object
  sh_ = new BDF2::SolutionHistory(2, map);
  state.init_solution_history(sh_);
  
  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&*paramList_,this);

}


Teuchos::RCP<Teuchos::ParameterList> BDF1Dae::getNonconstParameterList() {
  return paramList_;
}


Teuchos::RCP<Teuchos::ParameterList> BDF1Dae::unsetParameterList() {
  Teuchos::RCP<Teuchos::ParameterList> paramList = paramList_;
  paramList_ = Teuchos::null;
  return paramList;
}


Teuchos::RCP<const Teuchos::ParameterList> BDF1Dae::getParameterList() const {
  return paramList_;
}

Teuchos::RCP<const Teuchos::ParameterList> BDF1Dae::getValidParameters() const {
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;
  static RCP<const ParameterList> validParams;
  if (is_null(validParams)) {
    RCP<ParameterList>  pl = Teuchos::rcp(new ParameterList("steady time integrator"));
    Teuchos::setIntParameter("steady max iterations",
                             5,
                             "If during the steady state calculation, the number of iterations of the nonlinear solver exceeds this number, the subsequent time step is reduced.",
                             &*pl);
    Teuchos::setIntParameter("steady min iterations",
                             2,
                             "If during the steady state calculation, the number of iterations of the nonlinear solver exceeds this number, the subsequent time step is increased.",
                             &*pl);
    Teuchos::setIntParameter("steady limit iterations",
                             12,
                             "If during the steady state calculation, the number of iterations of the nonlinear solver exceeds this number, the current time step is reduced and the current time step is repeated.",
                             &*pl);
    Teuchos::setDoubleParameter("steady nonlinear tolerance",
                                0.1,
                                "The tolerance for the nonlinear solver during the steady state computation.",
                                &*pl);
    Teuchos::setDoubleParameter("steady time step reduction factor",
                                0.5,
                                "When time step reduction is necessary during the steady calculation, use this factor.",
                                &*pl);
    Teuchos::setDoubleParameter("steady time step increase factor",
                             1.2,
                             "When time step increase is possible during the steady calculation, use this factor.",
                             &*pl);
    Teuchos::setDoubleParameter("steady max time step",
                              1.0,
                              "The maximum allowed time step.",
                              &*pl);
    Teuchos::setIntParameter("steady max preconditioner lag iterations",
                              1,
                              "The maximum number of time steps that the preconditioner is allowed to be lagged.",
                              &*pl);
    Teuchos::setupVerboseObjectSublist(&*pl);
    validParams = pl;
  }
  return validParams;
}



void BDF1Dae::commit_solution(const double h, const Epetra_Vector& u) {

  double t = h + state.uhist->most_recent_time();

  state.uhist->record_solution(t, u);

  state.hlast = h;
  state.seq++;
  state.freeze_count = std::max<int>(0, state.freeze_count-1);

  state.hmin = std::min<double>(h, state.hmin);
  state.hmax = std::max<double>(h, state.hmax);

}


void BDF1Dae::select_step_size(const std::vector<double>& dt, const double perr, double& h) {
  const double tol = 0.001;

  ASSERT(dt.size() == 3);

  double a = 0.5*dt[0]*dt[1]*dt[2]/std::max<double>(perr,0.001);
  h = dt[0];

  double dh;
  do
  {
    double phi  = h*(h + dt[0])*(h + dt[1]) - a;
    double dphi = (2.0*h + dt[0])*(h + dt[1]) + h*(h + dt[0]);

    dh = phi / dphi;
    h = h - dh;
  }
  while (abs(dh) / h >= tol);

}




void BDF1Dae::set_initial_state(const double t, const Epetra_Vector& x, const Epetra_Vector& xdot) {
  // ASSERT(x.Map().PointSameAs( xdot.Map() ) );

  state.uhist->flush_history(t, x, xdot);
  state.seq = 0;
  state.usable_pc = false;
  state.pclagcount = 0;
}



void BDF1Dae::bdf1_step(double h, Epetra_Vector& u, double& hnext) {
  
  double tlast = state.uhist->most_recent_time();
  double tnew = tlast + h;
  
  bool fresh_pc = false;
  
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab
  
  hnext = h;

  // save the initial value so that we can restore the initial state 
  // in case we need to bail
  Epetra_Vector usav(map);
  usav = u;
  

  if (h < state.hmin) {
    std::cout << h << " " << state.hmin << std::endl;

    std::string msg = "BDF1 failed: Time step crash";
    Errors::Message m(msg);
    Exceptions::amanzi_throw(m);    
  }
  
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
    *out << "BDF1 step " << state.seq+1 << " T = "<< tlast << " H = " << h << std::endl;
  }
  
  // Predicted solution (initial value for the nonlinear solver)
  Epetra_Vector up(map);

  if ( state.uhist->history_size() > 1) {
    state.uhist->interpolate_solution(tnew,  up);
  } else  {
    up = u;
  }


  // u at the start of the time step
  Epetra_Vector u0(map);
  u0 = u;
  
  if (state.pclagcount > state.maxpclag) {
    state.usable_pc = false;
    state.pclagcount = 0;
  }
  
  if (!state.usable_pc) {
    // update the preconditioner (we use the predicted solution at the end of the time step)
    state.updpc_calls++;
    state.pclagcount++;
    int errc = 0;
    fn.update_precon (tnew, up, h, errc);
    if (errc != 0) {
      std::string msg = "BDF1 failed: error while updating the preconditioner.";
      Errors::Message m(msg);
      Exceptions::amanzi_throw(m);        
    }
    state.usable_pc = true;
  }

  //  Solve the nonlinear BCE system.
  u = up; // Initial solution guess is the predictor.
  try {
    solve_bce ( tnew, h, u0, u);
  }
  catch (int itr) { 
    // we end up in here either if the solver took too many iterations, 
    // or if it took too few
    if (itr > state.maxitr && itr <= state.mitr) { 
      hnext = std::min(state.hred * h, state.hlimit);
      state.usable_pc = false;
    } else if (itr < state.minitr) {
      hnext = std::min( state.hinc * h, state.hlimit);
      state.usable_pc = false;
    } else if (itr > state.mitr) {
      u = usav; // restore the original u
      state.usable_pc = false;
      hnext = std::min(h, state.hlimit);
      throw itr;
    } 
  }
}


void BDF1Dae::solve_bce(double t, double h, Epetra_Vector& u0, Epetra_Vector& u) {

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  fpa->nka_restart();

  int itr = 0;

  Epetra_Vector du(map);
  Epetra_Vector u_tmp(map);


  Teuchos::RCP<NOX::Epetra::Vector> preconditioned_f =
      Teuchos::rcp(new NOX::Epetra::Vector(du, NOX::ShapeCopy));

  do {
    double error;

    // Check for too many nonlinear iterations.
    if (itr >= state.mitr) {
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "AIN BCE solve failed " << itr << " iterations (max), error = " << error << std::endl;
      }
      throw itr;
    }

    itr++;

    // Evaluate the preconditioned nonlinear function.
    state.pcfun_calls++;

    // compute u_tmp = (u-u0)/h
    u_tmp = u;
    u_tmp.Update(-1.0/h,u0,1.0/h);

    fn.fun(t, u, u_tmp, du);

    fn.precon(du, u_tmp);

    // Accelerated correction.
    *preconditioned_f  = u_tmp;        // copy preconditioned functional into appropriate data type
    NOX::Epetra::Vector nev_du(du, NOX::ShapeCopy);  // create a vector for the solution

    fpa->nka_correction(nev_du, preconditioned_f);

    du = nev_du.getEpetraVector();  // copy result into an Epetra_Vector.

    // Next solution iterate and error estimate.
    // FORTRAN:  u  = u - du
    u.Update(-1.0,du,1.0);

    // Check the solution iterate for admissibility.
    if ( ! fn.is_admissible(u) ) { // iterate is bad; bail.
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "AIN BCE solve FAILED: inadmissible solution iterate: itr=" << itr << std::endl;
      }
      throw std::string("solution iterate is inadmissible"); 
    }

    error = fn.enorm(u, du);
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
      *out << itr << ": error = " << error << std::endl;
    }

    if (error > state.hlimit) {
      // the solver threatening to diverge
      throw state.mitr+1;
    }

    // Check for convergence
    if (error < state.ntol)   {
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "AIN BCE solve succeeded: " << itr << " iterations, error = "<< error << std::endl;
      }

      if ((itr < state.minitr) || (itr > state.maxitr)) {
        throw itr;
      }
      return;
    }
  }
  while (true);

}


void BDF1Dae::write_bdf1_stepping_statistics() {
  std::ostringstream oss;

  oss.flush();
  oss.setf(ios::scientific, ios::floatfield);

  oss << "STEP=";
  oss.fill('0');
  oss.width(5);
  oss << right << state.seq;
  oss << " T=";
  oss.precision(5);
  oss.width(11);
  oss << left << state.uhist->most_recent_time();
  oss << " H=";
  oss.precision(5);
  oss.width(11);
  oss << left << state.hlast;
  oss << " NFUN:NPC=";
  oss.fill('0');
  oss.width(4);
  oss << right << state.pcfun_calls;
  oss << ":";
  oss.fill('0');
  oss.width(4);
  oss << right << state.updpc_calls;
  oss << " NPCF:NNR:NNF:NSR=";
  oss.fill('0');
  oss.width(4);
  oss << right << state.updpc_failed;
  oss << ":";
  oss.fill('0');
  oss.width(2);
  oss << right << state.retried_bce;
  oss << ":";
  oss.fill('0');
  oss.width(2);
  oss << right << state.failed_bce;
  oss << ":";
  oss.fill('0');
  oss.width(2);
  oss << right << state.rejected_steps;


  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true))  {
    *out << oss.str() << std::endl;
  }
}

}
