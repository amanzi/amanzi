
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "NOX_Epetra_Vector.H"

#include "BDF2_TI.hh"
#include "SolutionHistory.hh"
#include "BDF_FnBase.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"


namespace BDF2 {

  const double Dae::RMIN = 0.25;
  const double Dae::RMAX = 4.0;
  const double Dae::MARGIN = 3.0;

  Dae::Dae(fnBase& fn_, const Epetra_BlockMap& map_) :
    fn(fn_), map(map_)
  {
    // set the line prefix for output
    this->setLinePrefix("BDF2::Dae           ");
    
    // make sure that the line prefix is printed
    this->getOStream()->setShowLinePrefix(true);
  }

  Dae::~Dae() {
    delete fpa;
    delete sh_;
  }


  void Dae::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
  {
    TEUCHOS_TEST_FOR_EXCEPT(is_null(paramList));

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
    state.mitr = paramList_->get<int>("Nonlinear solver max iterations");
    state.ntol = paramList_->get<double>("Nonlinear solver tolerance");
    
    mtries = paramList_->get<int>("Maximum number of BDF tries",20);

    int maxv = state.mitr-1;
    int mvec = paramList_->get<int>("NKA max vectors");
    maxv = std::min<int>(maxv,mvec);
    
    // Initialize the FPA structure.
    // first create a NOX::Epetra::Vector to initialize nka with
    NOX::Epetra::Vector init_vector( Epetra_Vector(map), NOX::ShapeCopy );
    double vtol = paramList_->get<double>("NKA drop tolerance");
    fpa = new nka(maxv, vtol, init_vector); 
    
    // create the solution history object
    sh_ = new SolutionHistory(3, map);
    state.init_solution_history(sh_);    

    // Read the sublist for verbosity settings.
    Teuchos::readVerboseObjectSublist(&*paramList_,this);
  }


  Teuchos::RCP<Teuchos::ParameterList> Dae::getNonconstParameterList()
  {
    return paramList_;
  }
  
  
  Teuchos::RCP<Teuchos::ParameterList> Dae::unsetParameterList()
  {
    Teuchos::RCP<Teuchos::ParameterList> paramList = paramList_;
    paramList_ = Teuchos::null;
    return paramList;
  }
  
  
  Teuchos::RCP<const Teuchos::ParameterList> Dae::getParameterList() const
  {
    return paramList_;
  }

  Teuchos::RCP<const Teuchos::ParameterList> Dae::getValidParameters() const
  {
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
  

  
  void Dae::commit_solution(const double h, const Epetra_Vector& u)
  {

    double t = h + state.uhist->most_recent_time();
    
    state.uhist->record_solution(t, u);
    
    state.hlast = h;
    state.seq++;
    state.freeze_count = std::max<int>(0, state.freeze_count-1);
    
    state.hmin = std::min<double>(h, state.hmin);
    state.hmax = std::max<double>(h, state.hmax);
  }


  void Dae::select_step_size(const std::vector<double>& dt, const double perr, double& h)
  {
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


  

  void Dae::set_initial_state(const double t, const Epetra_Vector& x, const Epetra_Vector& xdot)
  {
    // ASSERT(x.Map().PointSameAs( xdot.Map() ) );

    state.uhist->flush_history(t, x, xdot);
    state.seq = 0;  
  }

  
  void Dae::bdf2_step(double& h, double hmin, Epetra_Vector& u, double& hnext)
  {

    ASSERT(hmin>=0.0);
    ASSERT(hmin<=h);
    ASSERT(mtries>=1);

    using Teuchos::OSTab;
    Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab    

    int tries = 0;
    do
      {
	tries++;
	
	// Check for too many attempts at a single step.
	if (tries > mtries) 
	  {
	    Errors::Message m("BDF2 step failed");
	    Exceptions::amanzi_throw(m);
	  }
	
	// Check for a too-small step size.
	if (h < hmin) 
	  {
	    Errors::Message m("BDF2 step size too small");
	    Exceptions::amanzi_throw(m);
	  }	    

	// Attempt a BDF2 step.
	int errc = 0;
	bdf2_step_simple (h, u, hnext, errc, true);
	if (errc == 0) return;

	// Step failed; try again with the suggested step size.
	if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))
	  {
	    *out << "Changing H by a factor of " << hnext/h << std::endl;
	  }
	h = hnext;
      }
    while (true);

  }


  
  void Dae::bdf2_step_simple(double h, Epetra_Vector& u, double& hnext, int& errc, bool ectrl)
  {
    ASSERT(state.seq>=0);
    ASSERT(h>0);

    double mitr_sav;

    switch (state.seq) 
      {
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


  void Dae::trap_step_one(double h, Epetra_Vector& u, double& hnext, int& errc)
  {

    double tlast = state.uhist->most_recent_time();
    double t = tlast + h;
    double etah = 0.5 * h;
    double t0 = tlast + etah;

    using Teuchos::OSTab;
    Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab
    
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))
      {
	*out << "Trap step " << state.seq+1 << " tlast = " << tlast << " h = " << h << " etah = " << etah << std::endl;
      }
    
    
    // Predicted solution and base point for the BCE step.
    Epetra_Vector u0(map);

    state.uhist->interpolate_solution (t,  u, 1);
    state.uhist->interpolate_solution (t0, u0, 1);

    // Update the preconditioner.
    state.updpc_calls++;
    fn.update_precon (t, u, etah, errc); 
    if (errc != 0)
      {
	state.updpc_failed++;
	if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))
	  {
	    *out << "Preconditioner update failed at T = " << t << "  ETAH = " << etah << std::endl;
	  }
	hnext = 0.1 * h; // want to quickly find an acceptably small step size.
	errc = -1;
	return;
      }
    
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
      {
	*out << "Preconditioner updated at T = " << t << " etah = " << etah << std::endl;
      }
    state.hpc = etah;
    state.usable_pc = true;

    // Solve the nonlinear BCE system.
    solve_bce (t, etah, u0, u, errc);
    if (errc != 0) 
      {
	state.failed_bce++;
	hnext = 0.1 * h;  // want to quickly find an acceptably small step size.
	state.freeze_count = 1;
	errc = 1;
	return;
      }

    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
      {
	*out << "Step accepted: no local error control" << std::endl;
      }
    state.freeze_count = 2;
    hnext = h;
    errc = 0;
   
  }



  void Dae::bdf2_step_gen(double h, Epetra_Vector& u, double& hnext, int& errc, bool ectrl)
  {


    double tlast = state.uhist->most_recent_time();
    double t = tlast + h;
    double eta = (state.hlast + h) / (state.hlast + 2.0 * h);
    double etah = eta * h;
    double t0 = tlast + (1.0 - eta)*h;

    using Teuchos::OSTab;
    Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
      {
	*out << "BDF2 step " << state.seq+1 << " T = "<< tlast << " H = " << h << " ETAH = " << etah << std::endl;
      }
    bool fresh_pc = false;

    // If the PC step size is too different than the current step size we tag
    // it as unusable in order to preempt a possible nonlinear solve failure.
    if (state.usable_pc) 
      {
	if (state.hpc/etah > 1.0 + MARGIN) state.usable_pc = false;
	if (etah/state.hpc > 1.0 + MARGIN) state.usable_pc = false;
      }

    // Predicted solution and base point for BCE step.
    Epetra_Vector up(map);
    Epetra_Vector u0(map);
    
    state.uhist->interpolate_solution(t,  up, 2);
    state.uhist->interpolate_solution(t0, u0, 1);

    // check the predicted solution for admissibility 
    // and update preconditioner/halve time step if not
    if ( ! fn.is_admissible(up) )
      { 
	state.updpc_calls++;
	fn.update_precon(t, up, etah, errc);
	if ( errc != 0) // update failed
	  {
	    state.updpc_failed++;
	    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
	      {
		*out << "Preconditioner update FAILED at T=" << t << ", ETAH=" << etah << std::endl;
	      }
	    hnext = 0.5 * h;
	    
	    errc = -1;
	    return;
	  }
	if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
	  {
	    *out << "Preconditioner updated at T=" << t << ", ETAH=" << etah << std::endl;
	  }
	
        state.hpc = etah;
        state.usable_pc = true;
	fresh_pc = true;
      }

    // BCE loop:
    do
      {
	// Update the preconditioner if necessary.
	if (! state.usable_pc) 
	  {
	    state.updpc_calls++;
	    fn.update_precon (t, up, etah, errc); 
	    if (errc != 0)  // update failed; cut h and return error condition.
	      { 
		state.updpc_failed++;
		if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
		  {
		    *out << "Preconditioner update failed, T = "<< t << " etah = " << etah << std::endl;
		  }
		hnext = 0.5 * h;
		errc = -1;
		return;

	      }
	    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
	      {
		*out << "Preconditioner updated at T = "<< t << " etah = " << etah <<std::endl;
	      }
	    state.hpc = etah;
	    state.usable_pc = true;
	    fresh_pc = true;
	  }
	
	//  Solve the nonlinear BCE system.
	u = up; // Initial solution guess is the predictor.
	solve_bce ( t, etah, u0, u, errc );
	if (errc == 0) break;  // leave the do-while loop, the BCE step was successful.

	if (fresh_pc) // preconditioner was fresh; cut h and return error condition.
	  {
	    state.failed_bce++;
	    hnext = 0.5 * h;
	    state.freeze_count = 1;
	    errc = 1;
	    return;
	  }
	else // update the preconditioner and retry the nonlinear solve.
	  {
	    state.retried_bce++;
	    state.usable_pc = false;
	    continue; // cycle
	  }
      }
    while (true);
	     

    //bool predictor_error = true;
    bool predictor_error = ectrl;

    if (predictor_error) 
      {
	
	// Predictor error control.
	// FORTRAN: u0 = u - up
	u0 = u;
	u0.Update(-1.0, up, 1.0);
	
	
	double perr = fn.enorm(u, u0);
	if (true) //  (perr < 4.0) // accept the step. 
	  { 
	    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
	      {

		*out << "Step accepted, perr = " << perr <<std::endl;
	      }
	    errc = 0;
	  }
	else  // reject the step; cut h and return error condition.
	  {
	    
	    state.rejected_steps++;
	    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
	      {
		*out << "Step REJECTED, perr = " << perr << std::endl;
	      }
	    hnext = 0.5 * h;
	    state.freeze_count = 1;
	    errc = 2;
	    return;
	  }

	// Select the next step size based on the predictor error and past step
	// size history, but don't change the step size by too great a factor.
	// FORTRAN:  dt(1) = h
	// FORTRAN:  dt(2:) = h + time_deltas(this%uhist)

	std::vector<double> dt(3);
	std::vector<double> dt0(2); // to store the time deltas from the solution history
	state.uhist->time_deltas(dt0);	
	
	dt[0] = h;
	dt[1] = h + dt0[0];
	dt[2] = h + dt0[1];

	select_step_size (dt, perr, hnext);
	
	hnext = std::max<double>(RMIN*h, std::min<double>(RMAX*h, hnext));
	if (state.freeze_count != 0) hnext = std::min<double>(h, hnext);
      }
    else
      {
	if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
	  {
	    *out << "Step accepted: no local error control" << std::endl;
	  }
	hnext = h;
	errc = 0;

      }
  }


  void Dae::solve_bce(double t, double h, Epetra_Vector& u0, Epetra_Vector& u, int& errc)
  {
    
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
    do
      {
	// Check for too many nonlinear iterations.
	if (itr >= state.mitr) 
	  {
	    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
	      {
		*out << "AIN BCE solve failed " << itr << " iterations (max), error = " << error << std::endl;
	      }
	    errc = 1;
	    return;
	  }
	
	itr++;

	// Evaluate the preconditioned nonlinear function.
	state.pcfun_calls++;
	
	// compute u_tmp = (u-u0)/h
	u_tmp = u;
	u_tmp.Update(-1.0/h,u0,1.0/h);
	
	fn.fun(t, u, u_tmp, du, h);

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
	if ( ! fn.is_admissible(u) ) // iterate is bad; bail.
	  {
	    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
	      {
		*out << "AIN BCE solve FAILED: inadmissible solution iterate: itr=" << itr << std::endl;
		errc = 2;
		return;
	      }
	  }
	
	error = fn.enorm(u, du); 
	if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
	  {
	    *out << itr << ": error = " << error << std::endl;
	  }

	// Check for convergence.
	if (((error < state.ntol) && (itr > 1)) || (error < 0.01 * state.ntol))
	  {
	    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH,true))	  
	      {
		*out << "AIN BCE solve succeeded: " << itr << " iterations, error = "<< error << std::endl;
	      }
	    errc = 0;
	    return;
	  }
      }
    while (true);

  }


  void Dae::write_bdf2_stepping_statistics()
  {
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
    
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true))	
      {
	*out << oss.str() << std::endl;
      }
  }

}
