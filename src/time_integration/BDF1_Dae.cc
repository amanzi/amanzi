
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

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

#include "Timer.hh"


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

  state.maxpclag = paramList_->get<int>("max preconditioner lag iterations");


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
    RCP<ParameterList>  pl = Teuchos::rcp(new ParameterList("time integrator"));
    Teuchos::setIntParameter("max iterations",
                             5,
                             "If during the steady state calculation, the number of iterations of the nonlinear solver exceeds this number, the subsequent time step is reduced.",
                             &*pl);
    Teuchos::setIntParameter("min iterations",
                             2,
                             "If during the steady state calculation, the number of iterations of the nonlinear solver exceeds this number, the subsequent time step is increased.",
                             &*pl);
    Teuchos::setIntParameter("limit iterations",
                             12,
                             "If during the steady state calculation, the number of iterations of the nonlinear solver exceeds this number, the current time step is reduced and the current time step is repeated.",
                             &*pl);
    Teuchos::setDoubleParameter("nonlinear tolerance",
                                0.1,
                                "The tolerance for the nonlinear solver during the steady state computation.",
                                &*pl);
    Teuchos::setDoubleParameter("time step reduction factor",
                                0.5,
                                "When time step reduction is necessary during the steady calculation, use this factor.",
                                &*pl);
    Teuchos::setDoubleParameter("time step increase factor",
                             1.2,
                             "When time step increase is possible during the steady calculation, use this factor.",
                             &*pl);
    Teuchos::setDoubleParameter("max time step",
                              1.0,
                              "The maximum allowed time step.",
                              &*pl);
    Teuchos::setIntParameter("max preconditioner lag iterations",
                              1,
                              "The maximum number of time steps that the preconditioner is allowed to be lagged.",
                              &*pl);
    Teuchos::setDoubleParameter("error abs tol",
                              1.0,
                              "Absolute error prefactor.",
                              &*pl);
    Teuchos::setDoubleParameter("error rel tol",
                              1.0,
                              "Relative error prefactor.",
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
  

  if (h < state.hmin) {
    *out << h << " " << state.hmin << std::endl;

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

//     solve_bce(tnew, h, u0, u);
    solve_bce_jfnk(tnew, h, u0, u);
//     exit(0);
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
}


void BDF1Dae::solve_bce(double t, double h, Epetra_Vector& u0, Epetra_Vector& u) {

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab
  
  Timer ttotal;
  
  ttotal.start();

  fpa->nka_restart();

  int itr = 0;

  Epetra_Vector du(map);
  Epetra_Vector u_tmp(map);


  Teuchos::RCP<NOX::Epetra::Vector> preconditioned_f =
      Teuchos::rcp(new NOX::Epetra::Vector(du, NOX::ShapeCopy));

  double error(0.0);
  double du_norm(0.0), previous_du_norm(0.0);
  int divergence_count(0);

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
    if (itr%(state.maxpclag+1)==0) {
      Timer t1;
      t1.start();
      fn.update_precon (t, u, h, errc);
      t1.stop();
      std::cout << "Preconditioner, if necessary: " << t1 << std::endl;
    }
 

    // iteration counter
    itr++;

    // count the number of nonlinear function evaluations
    state.pcfun_calls++;

    // compute u_tmp = (u-u0)/h
    u_tmp = u;
    u_tmp.Update(-1.0/h, u0, 1.0/h);

    // evaluate nonlinear functional
//     DEBUG: Nathan
    Timer t1;
    t1.start();
    fn.fun(t, u, u_tmp, du, h);
    t1.stop();
    std::cout << "evaluate nonlinear functional: " << t1 << std::endl;
    
    // apply preconditioner to the nonlinear residual
    // DEBUG: Nathan
    Timer t2;
    t2.start();
    fn.precon(du, u_tmp);
    t2.stop();
    std::cout << "apply preconditioner to the nonlinear residual: " << t2 << std::endl;   
    
    // stuff the preconditioned residual into a NOX::Epetra::Vector
    *preconditioned_f = u_tmp;  // copy preconditioned functional into appropriate data type
    NOX::Epetra::Vector nev_du(du, NOX::ShapeCopy);  // create a vector for the solution
    
    

    // compute the accellerated correction
    fpa->nka_correction(nev_du, preconditioned_f);
    
    

    // copy result into an Epetra_Vector.
    du = nev_du.getEpetraVector();  


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
      throw state.mitr+1;
    }

    // protect against divergencve of the nonlinear solver
    if (itr > 1 && du_norm >= previous_du_norm) {
      // the solver threatening to diverge
      ++divergence_count;

      // if it does not recover quickly, abort
      if (divergence_count == 3) {
	*out << "Nonlinear solver is starting to diverge, cutting current time step." << std::endl;
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
        *out << "AIN BCE solve succeeded: " << itr << " iterations, error = "<< error << std::endl;

	ttotal.stop();
  
        std::cout << "nonlinear solver takes: " << ttotal << std::endl;
	
	//    // Print solution
	char file_name[25];
	FILE *ifp;
	int MyPID = 0;
	
	(void) sprintf(file_name, "output_nka.%d",MyPID);
	ifp = fopen(file_name, "w");
	for (int i=0; i<u.MyLength(); i++)
	fprintf(ifp, "%d %E\n", i, u[i]);
	fclose(ifp);

	exit(0);

      }

      if ((itr < state.minitr) || (itr > state.maxitr)) {
        throw itr;
      }

      if (divergence_count > 0) throw state.maxitr+1;
      
      
  
      return;
    }
  }
  while (true);
  
  

}


void BDF1Dae::solve_bce_jfnk(double t, double h, Epetra_Vector& u0, Epetra_Vector& u) {
	
  Timer ttotal;
  
  ttotal.start();
	
  int NumMyElements = u.MyLength();
//   for (int i=0;i<NumMyElements;i++) u[i] = 1;
	
	
  NOX::Epetra::Vector nox_u(u);
  
//   nox_u =  u;
  	
  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *nlParamsPtr.get();

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");
  //nlParams.set("Nonlinear Solver", "Trust Region Based");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
//   printParams.set("MyPID", MyPID); 
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
//   printParams.set("Output Information", 
// 			NOX::Utils::OuterIteration + 
// 			NOX::Utils::OuterIterationStatusTest + 
// 			NOX::Utils::InnerIteration +
// 			NOX::Utils::Parameters + 
// 			NOX::Utils::Details + 
// 			NOX::Utils::Warning);

  // Create printing utilities
  NOX::Utils utils(printParams);

  // Sublist for line search 
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "Full Step");
  //searchParams.set("Method", "Interval Halving");
  //searchParams.set("Method", "Polynomial");
  //searchParams.set("Method", "NonlinearCG");
  //searchParams.set("Method", "Quadratic");
  //searchParams.set("Method", "More'-Thuente");

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
//  dirParams.set("Method", "Modified-Newton");
//  Teuchos::ParameterList& newtonParams = dirParams.sublist("Modified-Newton");
//    newtonParams.set("Max Age of Jacobian", 2);
  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");
//     newtonParams.set("Forcing Term Method", "Type 1");
    //newtonParams.set("Forcing Term Method", "Type 2");
//     newtonParams.set("Forcing Term Minimum Tolerance", 1.0e-4);
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
  lsParams.set("Max Iterations", 800);  
  lsParams.set("Tolerance", 1e-4); 
  lsParams.set("Preconditioner", "None");
//   lsParams.set("Preconditioner", "Ifpack");
  lsParams.set("Max Age Of Prec", 5); 

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  // NOX_Epetra_Interface
  
//   AmanziFlow::Flow_PK* flow_problem = dynamic_cast<AmanziFlow::Flow_PK*> &fn;
  
  
  
  const Teuchos::RCP<NOX::Epetra::Interface::Required> interface = 
    Teuchos::rcp(new AmanziFlow::Interface_NOX(&fn, u0, t, h));

  // Create the Epetra_RowMatrix.  Uncomment one or more of the following:
  // 1. User supplied (Epetra_RowMatrix)
  //Teuchos::RCP<Epetra_RowMatrix> Analytic = Problem.getJacobian();
  // 2. Matrix-Free (Epetra_Operator)
  Teuchos::RCP<NOX::Epetra::MatrixFree> MF = 
    Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface, nox_u));
//   3. Finite Difference (Epetra_RowMatrix)
  Teuchos::RCP<NOX::Epetra::FiniteDifference> FD = 
      Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams, interface, nox_u));

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = MF;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = FD;
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
						      iJac,  MF, iPrec, FD, nox_u));

  // Create the Group
  Teuchos::RCP<NOX::Epetra::Group> grp =
    Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, nox_u, 
					linSys)); 
//        
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormF> relresid = 
    Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), 1.0e-7));
  Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
    Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  converged->addStatusTest(absresid);
  converged->addStatusTest(relresid);
  converged->addStatusTest(wrms);
  converged->addStatusTest(update);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(200));
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
  NOX::StatusTest::StatusType status = solver->solve();
  
//   interface -> printTime();

  if (status == NOX::StatusTest::Converged){
    utils.out() << "Test Passed!" << endl;
    throw solver->getNumIterations();
  }
  else {
    utils.out() << "Nonlinear solver failed to converge!" << endl;
    throw state.maxitr+1;
   }

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
//   const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
    u = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
// 
//   // End Nonlinear Solver **************************************
    
//      // Output the parameter list
//   if (utils.isPrintType(NOX::Utils::Parameters)) {
//     utils.out() << endl << "Final Parameters" << endl
// 	 << "****************" << endl;
//     solver->getList().print(utils.out());
//     utils.out() << endl;
//   }
// // 
    
    
    
    
    
    ttotal.stop();
  
    std::cout << "nonlinear solver takes: " << ttotal << std::endl;
    
    
//    // Print solution
//   char file_name[25];
//   FILE *ifp;
//   int MyPID = 0;
//   
//   (void) sprintf(file_name, "output_jfnk.%d",MyPID);
//   ifp = fopen(file_name, "w");
//   for (int i=0; i<NumMyElements; i++)
//     fprintf(ifp, "%d %E\n", i, u[i]);
//   fclose(ifp);
// 
  exit(0);
	
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
