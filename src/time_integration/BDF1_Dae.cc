
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "InputParserIS-defaults.hh"

#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_Epetra_Vector.H"

#include "Interface_NOX.hpp"

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

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(paramList),std::runtime_error,"Parameter list is null");

  // Validate and set the parameter defaults.  Here, the parameters are
  // validated and the state of *this is not changed unless the parameter
  // validation succeeds.  Also, any validators that are defined for various
  // parameters are passed along so theu they can be used in extracting
  // values!


  //paramList->validateParametersAndSetDefaults(*this->getValidParameters(),0);
  paramList_ = paramList;

  // make sure that the parameter list is actually valid (this is probably redundant)
  //paramList_->validateParameters(*this->getValidParameters());

  // read the parameter list and initialize the class

  state.mitr = paramList_->get<int>("limit iterations");
  state.maxitr = paramList_->get<int>("max iterations"); 
  state.minitr = paramList_->get<int>("min iterations");
  state.ntol = paramList_->get<double>("nonlinear tolerance");
  state.hred = paramList_->get<double>("time step reduction factor");
  state.hinc = paramList_->get<double>("time step increase factor");
  state.hlimit = paramList_->get<double>("max time step");
  state.atol = paramList_->get<double>("error abs tol");
  state.rtol = paramList_->get<double>("error rel tol");
  state.damp = paramList_->get<double>("nonlinear iteration damping factor",1.0);
  state.uhist_size = paramList_->get<int>("nonlinear iteration initial guess extrapolation order",1);
  state.uhist_size++;

  state.maxpclag = paramList_->get<int>("max preconditioner lag iterations");
  state.currentpclag = state.maxpclag;

  state.nonlinear_solver = BDFNKA;
  
  max_divergence_count_ = paramList_->get<int>("max divergent iterations", MAX_DIVERGENT_ITERATIONS);

  std::string nstype = paramList_->get<std::string>("nonlinear solver", "NKA");
  if (nstype == "NKA") {
      state.nonlinear_solver = BDFNKA;
  } else if (nstype == "JFNK") {
      state.nonlinear_solver = BDFJFNK;
  } else {
       Errors::Message m("(native spec) bdf1 nonlinear solver must either be NKA or JFNK");
       Exceptions::amanzi_throw(m);
  }


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
  sh_ = new BDF2::SolutionHistory(state.uhist_size, map);
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
    RCP<ParameterList>  pl = Teuchos::rcp(new ParameterList("time integrator"));
    setIntParameter("max iterations",
                             5,
                             "If during the steady state calculation, the number of iterations of the nonlinear solver exceeds this number, the subsequent time step is reduced.",
                             &*pl);
    setIntParameter("min iterations",
                             2,
                             "If during the steady state calculation, the number of iterations of the nonlinear solver exceeds this number, the subsequent time step is increased.",
                             &*pl);
    setIntParameter("limit iterations",
                             12,
                             "If during the steady state calculation, the number of iterations of the nonlinear solver exceeds this number, the current time step is reduced and the current time step is repeated.",
                             &*pl);
    setDoubleParameter("nonlinear tolerance",
                                0.1,
                                "The tolerance for the nonlinear solver during the steady state computation.",
                                &*pl);
    setDoubleParameter("time step reduction factor",
                                0.5,
                                "When time step reduction is necessary during the steady calculation, use this factor.",
                                &*pl);
    setDoubleParameter("time step increase factor",
                             1.2,
                             "When time step increase is possible during the steady calculation, use this factor.",
                             &*pl);
    setDoubleParameter("max time step",
                              1.0,
                              "The maximum allowed time step.",
                              &*pl);
    setIntParameter("max preconditioner lag iterations",
                              1,
                              "The maximum number of time steps that the preconditioner is allowed to be lagged.",
                              &*pl);
    setDoubleParameter("error abs tol",
                              1.0,
                              "Absolute error prefactor.",
                              &*pl);
    setDoubleParameter("error rel tol",
                              1.0,
                              "Relative error prefactor.",
                              &*pl);
    setIntParameter("max divergent iterations",
		              3,
		    "Maximum divergent nonlinear iterations the bdf1 time iterator will tolerate",
		              &*pl);
    setDoubleParameter("nka damping factor",
		       1.0,
		       "NKA damping factor.",
		       &*pl);
    setupVerboseObjectSublist(&*pl);
    validParams = pl;
  }
  return validParams;
}



void BDF1Dae::commit_solution(const double h, const Epetra_Vector& u) {

  double t = h + state.uhist->most_recent_time();

  state.uhist->record_solution(t, u);

  state.hlast = h;
  state.seq++;
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

  fn.update_norm(state.rtol,state.atol);
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
  
  if (h < 10.0 * state.eps * tlast) {
    std::string msg = "BDF1 failed: Time step crash";
    Errors::Message m(msg);
    Exceptions::amanzi_throw(m);    
  }
  
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
    *out << "BDF1 step " << state.seq+1 << " T = "<< tlast << " H = " << h << std::endl;
  }
  
  // Predicted solution (initial value for the nonlinear solver)
  Epetra_Vector up(map);

  if (state.uhist->history_size() > 1) {
    state.uhist->interpolate_solution(tnew,  up);
  } else {
    up = u;
  }

  // allow the pk to modify the predictor
  bool predictor_modified = fn.modify_predictor(h,up);

  // u at the start of the time step
  Epetra_Vector u0(map);
  u0 = u;
  
  if (state.pclagcount > state.maxpclag) {
    state.usable_pc = false;
    state.pclagcount = 0;
  }
  
  //  Solve the nonlinear BCE system.
  u = up;  // Initial solution guess is the predictor.
  
  
    
  try {
      if (state.nonlinear_solver == BDFNKA) {
            solve_bce(tnew, h, u0, u);            
      }
      else if (state.nonlinear_solver == BDFJFNK) {
             solve_bce_jfnk(tnew, h, u0, u);
      }     
  }
  catch (int itr) { 
    // we end up in here either if the solver took too many iterations, 
    // or if it took too few        
    if (itr > state.maxitr && itr <= state.mitr) { 
      hnext = std::min(state.hred * h, state.hlimit);
      state.usable_pc = false;
    } else if (itr < state.minitr) {
      hnext = std::min(state.hinc * h, state.hlimit);
      state.usable_pc = false;
    } else if (itr > state.mitr) {
      u = usav;  // restore the original u
      state.usable_pc = false;
      hnext = std::min(h, state.hlimit);
      throw itr;
    } 
  }
  
//   utmp.Update(-1, u, 1);
//   
//   double val;
//   
//     std::cout<<"After solve_bce??? "<< val<<"\n";
//   exit(0);
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

  double error(0.0);
  double du_norm(0.0), previous_du_norm(0.0);
  int divergence_count(0);

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
    *out << "BDF1: preconditioner lag is " << state.currentpclag  <<std::endl;
  }  

  do {
    // Check for too many nonlinear iterations.
    if (itr > state.mitr) {
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "AIN BCE solve failed " << itr << " iterations (max), error = " << error << std::endl;
      }
      throw itr;
    }

    // update the preconditioner if necessary
    // DEBUG: Nathan
    int errc(0);
    if (itr%(state.currentpclag+1)==0) {
      fn.update_precon (t, u, h, errc);
    }
 

    // iteration counter
    itr++;

    // count the number of nonlinear function evaluations
    state.pcfun_calls++;

    // compute u_tmp = (u-u0)/h
    u_tmp = u;
    u_tmp.Update(-1.0/h, u0, 1.0/h);

    // evaluate nonlinear functional
    fn.fun(t, u, u_tmp, du, h);

    // apply preconditioner to the nonlinear residual
    fn.precon(du, u_tmp);
    
    // stuff the preconditioned residual into a NOX::Epetra::Vector
    *preconditioned_f = u_tmp;  // copy preconditioned functional into appropriate data type
    NOX::Epetra::Vector nev_du(du, NOX::ShapeCopy);  // create a vector for the solution

    // compute the accellerated correction
    if (itr>1) {
      fpa->nka_correction(nev_du, preconditioned_f, state.damp);
    } else {
      fpa->nka_correction(nev_du, preconditioned_f);
    }

    // copy result into an Epetra_Vector.
    du = nev_du.getEpetraVector();  
    
    // apply damping
    du.Scale(state.damp);

    // Check the solution iterate for admissibility.
    if ( ! fn.is_admissible(u) ) { // iterate is bad; bail.
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "AIN BCE solve FAILED: inadmissible solution iterate: itr=" << itr << std::endl;
      }
      throw std::string("solution iterate is inadmissible"); 
    }

    // make sure that we do not diverge and cause numerical overflow
    // we use inf norms here
    previous_du_norm = du_norm;
    du.NormInf(&du_norm);

    // protect against floating point overflow
    if (itr > 1 && du_norm > 1000.0 * previous_du_norm) {
      *out << itr << ": error (infinity norm) " << du_norm << std::endl;
      *out << "Nonlinear solver is threatening to overflow, cutting current time step." << std::endl;
      state.currentpclag = 0;
      throw state.mitr+1;
    }

    // protect against divergencve of the nonlinear solver
    if (itr > 1 && du_norm >= previous_du_norm) {
      // the solver threatening to diverge
      ++divergence_count;

      // if it does not recover quickly, abort
      if (divergence_count == max_divergence_count_) {
	*out << "Nonlinear solver is starting to diverge, cutting current time step." << std::endl;
	state.currentpclag = 0;
	throw state.mitr+1;
      }
    } else {
      divergence_count = 0;
    }

    // next solution iterate and error estimate.
    //   u  = u - du
    u.Update(-1.0,du,1.0);

    // compute error norm for the purpose of convergence testing, we use the
    // norm provided in the model evaluator
    error = fn.enorm(u, du);
    
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
      *out << itr << ": error = " << error << std::endl;
    }

    // Check for convergence
    if (error < state.ntol)   {
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        *out << "AIN BCE solve succeeded: " << itr << " iterations, error = "<< error <<std::endl;
      }

      if (itr < state.minitr) {
	state.currentpclag = std::min(state.currentpclag+1, state.maxpclag);
	throw itr;
      }
      if (itr > state.maxitr) {
	state.currentpclag = std::max(state.currentpclag-1, 0);
        throw itr;
      }

      return;
    }
  }
  while (true);
}


void BDF1Dae::solve_bce_jfnk(double t, double h, Epetra_Vector& u0, Epetra_Vector& u) {

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

//   u0.Print(std::cout);
//   u.Print(std::cout);
  
  NOX::Epetra::Vector nox_u(u);

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *nlParamsPtr.get();

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");
  // nlParams.set("Nonlinear Solver", "Trust Region Based");

  nlParams.set("Jacobian Operator", "Matrix-Free");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) {
        //   printParams.set("MyPID", MyPID);

//         printParams.set("Linear Solver Details", true);
        printParams.set("Output Information",
             NOX::Utils::OuterIteration + 
             NOX::Utils::OuterIterationStatusTest + 
             NOX::Utils::InnerIteration +
//              NOX::Utils::Parameters + 
//              NOX::Utils::Details + 
             NOX::Utils::Warning);
  }
  else if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true)){
           printParams.set("Output Information",
                NOX::Utils::OuterIteration + 
//                 NOX::Utils::OuterIterationStatusTest + 
                NOX::Utils::InnerIteration +
                NOX::Utils::Warning);    
  }
  else if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)){
           printParams.set("Output Information",
             NOX::Utils::Warning);    
  }
  // Create printing utilities
  NOX::Utils utils(printParams);

  // Sublist for line search 
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "Full Step");
  //searchParams.set("Method", "Backtrack");
  //searchParams.set("Method", "Interval Halving");
  //searchParams.set("Method", "Polynomial");
  //searchParams.set("Method", "NonlinearCG");
  //searchParams.set("Method", "Quadratic");
  //searchParams.set("Method", "More'-Thuente");

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  //dirParams.set("Method", "Modified-Newton");
  //Teuchos::ParameterList& newtonParams = dirParams.sublist("Modified-Newton");
  //newtonParams.set("Max Age of Jacobian", 2);
  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  newtonParams.set("Forcing Term Method", "Constant");
  //newtonParams.set("Forcing Term Method", "Type 1");
  //newtonParams.set("Forcing Term Method", "Type 2");
  //newtonParams.set("Forcing Term Minimum Tolerance", 1.0e-4);
  //newtonParams.set("Forcing Term Maximum Tolerance", 0.1);
  //dirParams.set("Method", "Steepest Descent");
  //Teuchos::ParameterList& sdParams = dirParams.sublist("Steepest Descent");
    //sdParams.set("Scaling Type", "None");
    //sdParams.set("Scaling Type", "2-Norm");
    //sdParams.set("Scaling Type", "Quadratic Model Min");
    //dirParams.set("Method", "NonlinearCG");
  //Teuchos::ParameterList& nlcgParams = dirParams.sublist("Nonlinear CG");
    //nlcgParams.set("Restart Frequency", 2000);
    //nlcgParams.set("Precondition", "On");
    //nlcgParams.set("Orthogonalize", "Polak-Ribiere");
    //nlcgParams.set("Orthogonalize", "Fletcher-Reeves");

  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");  
  lsParams.set("Max Iterations", 5);  
  lsParams.set("Tolerance", 1e-2); 
  //lsParams.set("Preconditioner", "None");
  //lsParams.set("Preconditioner", "Ifpack");
  lsParams.set("Preconditioner", "User Defined");
  lsParams.set("Max Age Of Prec", 1); 

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  // NOX_Epetra_Interface
  

  const Teuchos::RCP<AmanziFlow::Interface_NOX> interface = 
    Teuchos::rcp(new AmanziFlow::Interface_NOX(&fn, &state, u0, t, h));
    
    interface->resetPrecLagCounter();
    interface->setPrecLag(5);

  // Create the Epetra_RowMatrix.  Uncomment one or more of the following:
  // 1. User supplied (Epetra_RowMatrix)
  //Teuchos::RCP<Epetra_RowMatrix> Analytic = Problem.getJacobian();
  // 2. Matrix-Free (Epetra_Operator)
  Teuchos::RCP<NOX::Epetra::MatrixFree> MF = 
    Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface, nox_u));
  // 3. Finite Difference (Epetra_RowMatrix)
  //Teuchos::RCP<NOX::Epetra::FiniteDifference> FD = 
  //    Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams, interface, nox_u));
    
///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

   Teuchos::RCP<AmanziFlow::JFNK_Preconditioner> Prec_Test = 
     Teuchos::rcp(new AmanziFlow::JFNK_Preconditioner(&fn));
///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = MF;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;

  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
                        iJac,  MF, iPrec, Prec_Test, nox_u));
//                        iReq, iJac, MF, nox_u));

  // Create the Group
  Teuchos::RCP<NOX::Epetra::Group> grp =
    Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, nox_u, 
                            linSys)); 
//
//   Create the convergence tests
//   Teuchos::RCP<NOX::StatusTest::NormF> absresid = 
//     Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-2));
//   Teuchos::RCP<NOX::StatusTest::NormF> converged = 
//     Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), 1.0e-8));
//   Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
//     Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-0));
//   Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
//     Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
//   Teuchos::RCP<NOX::StatusTest::Combo> converged =
//     Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
//   converged->addStatusTest(absresid);
//   converged->addStatusTest(relresid);
//   converged->addStatusTest(wrms);
//   converged->addStatusTest(update);
    
    Teuchos::RCP<AmanziFlow::PK_enorm> converged =
        Teuchos::rcp(new AmanziFlow::PK_enorm(&fn, state.ntol));
    
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(state.mitr));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);


  Teuchos::RCP<Teuchos::ParameterList> finalParamsPtr = nlParamsPtr;


  // Create the method
  Teuchos::RCP<NOX::Solver::Generic> solver = 
    NOX::Solver::buildSolver(grp, combo, finalParamsPtr);

    NOX::StatusTest::StatusType status;


    try {
        status = solver->solve();
    }
    catch (int ierr){
        if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true)) {
                *out <<"Error in JFNK solver "<<ierr<<std::endl;
                *out <<status<<std::endl;
        }
     }

//   interface -> printTime();


  if (status == NOX::StatusTest::Converged) {
    // if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) 
    utils.out() << "Test Passed!" << endl;
    
     // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
//   const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
    u = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

    write_bdf1_stepping_statistics();

    throw solver->getNumIterations();
  }
  else 
  {
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true)) utils.out() << "Nonlinear solver failed to converge!" << endl;
    throw state.mitr+1;
   }

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

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true))  {
    *out << oss.str() << std::endl;
  }
}

}
