#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputParserIS.hh"
#include "InputParserIS_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

/* ******************************************************************
* This routine has to be called after routine CreateStateList so that
* constant density will be properly initialized.
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateFlowList_(Teuchos::ParameterList* plist, int time_regime) 
{
  Errors::Message msg;
  Teuchos::OSTab tab = vo_->getOSTab();
  Teuchos::ParameterList flw_list;
  Teuchos::ParameterList *flow_list;

  if (plist->isSublist("Execution Control")) {
    Teuchos::ParameterList& exe_list = plist->sublist("Execution Control");

    if (exe_list.isParameter("Flow Model")) {
      std::string flow_model = exe_list.get<std::string>("Flow Model");

      // set up default values for some expert parameters
      std::string disc_method("MFD: Optimized for Sparsity");
      std::string prec_method("Linearized Operator");
      std::string rel_perm("Upwind: Amanzi");
      std::string update_upwind("every timestep");
      double atm_pres(ATMOSPHERIC_PRESSURE);
      std::string nonlinear_solver("NKA");
      bool modify_correction(false);

      // process expert parameters
      if (exe_list.isSublist("Numerical Control Parameters")) {
        Teuchos::ParameterList& ncp_list = exe_list.sublist("Numerical Control Parameters");

        if (ncp_list.isSublist("Unstructured Algorithm")) {
          Teuchos::ParameterList& ua_list = ncp_list.sublist("Unstructured Algorithm");
          if (ua_list.isSublist("Flow Process Kernel")) {
            Teuchos::ParameterList fl_exp_params = ua_list.sublist("Flow Process Kernel");
            if (fl_exp_params.isParameter("Discretization Method")) {
              disc_method = fl_exp_params.get<std::string>("Discretization Method");
            }
            if (fl_exp_params.isParameter("Relative Permeability")) {
              rel_perm = fl_exp_params.get<std::string>("Relative Permeability");
            }
            if (fl_exp_params.isParameter("atmospheric pressure")) {
              atm_pres = fl_exp_params.get<double>("atmospheric pressure");
            }
            if (fl_exp_params.isParameter("Preconditioning Strategy")) {
              prec_method = fl_exp_params.get<std::string>("Preconditioning Strategy");
            }
            if (fl_exp_params.isParameter("Use Picard")) {
              use_picard_ = fl_exp_params.get<bool>("Use Picard");
            }
          }
          if (ua_list.isSublist("Nonlinear Solver")) {
            Teuchos::ParameterList fl_exp_params = ua_list.sublist("Nonlinear Solver");
            if (fl_exp_params.isParameter("Nonlinear Solver Type")) {
              nonlinear_solver = fl_exp_params.get<std::string>("Nonlinear Solver Type");
            }
	    update_upwind = fl_exp_params.get<std::string>("update upwind frequency", "every timestep");
            modify_correction = fl_exp_params.get<bool>("modify correction", false);
          }
        }
      }

      // Newton method requires to overwrite other parameters.
      if (nonlinear_solver == std::string("Newton")) {
        disc_method = std::string("FV: Default");
        rel_perm = std::string("Upwind: Darcy Velocity");
	update_upwind = std::string("every nonlinear iteration");	
        if (vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
          *vo_->os() << vo_->color("yellow") << "FV enforces: \"Upwind: Darcy Velocity\" "
                     << ", \"modify correction\", and \"FV: Default\"." 
                     << vo_->reset() << std::endl;
        }
      }
      
      if (flow_model == "Single Phase" || flow_model == "Richards") {
        if (flow_model == "Single Phase") {
          Teuchos::ParameterList& darcy_problem = flw_list.sublist("Darcy problem");
          darcy_problem.sublist("VerboseObject") = CreateVerbosityList_(verbosity_level);
          darcy_problem.set<double>("atmospheric pressure", atm_pres);

          flow_list = &darcy_problem; // we use this below to insert sublists that are shared by Richards and Darcy
          flow_single_phase = true;
        }
        else if (flow_model == "Richards") {
          Teuchos::ParameterList& richards_problem = flw_list.sublist("Richards problem");

          Teuchos::ParameterList& upw_list = richards_problem.sublist("upwind");
          upw_list.set<std::string>("relative permeability", boost::to_lower_copy(rel_perm));
          upw_list.set<std::string>("upwind update", update_upwind);
          // "standard" is the most robust upwind method for variety of subsurface
          // scenarios. Note that "Upwind: Amanzi" requires "upwind method"="divk" 
          // to reproduce the same behavior on orthogonal meshes. 
          if (rel_perm == "Upwind: Amanzi") {
            upw_list.set<std::string>("upwind method", "divk");
            upw_list.sublist("upwind divk parameters").set<double>("tolerance", 1e-12);
          } else {
            upw_list.set<std::string>("upwind method", "standard");
            upw_list.sublist("upwind standard parameters").set<double>("tolerance", 1e-12);
          }
          // this one should come from the input file...
          richards_problem.sublist("VerboseObject") = CreateVerbosityList_(verbosity_level);
          richards_problem.set<double>("atmospheric pressure", atm_pres);

          // see if we need to generate a Picard list
          flow_list = &richards_problem; // we use this below to insert sublists that are shared by Richards and Darcy

          // insert the water retention models sublist (these are only relevant for Richards)
          Teuchos::ParameterList& water_retention_models = richards_problem.sublist("water retention models");
          water_retention_models = CreateWRM_List_(plist);
        }

        // insert operator sublist
        Teuchos::ParameterList op_list;
        op_list = CreateFlowOperatorList_(disc_method, prec_method, nonlinear_solver, rel_perm);
        flow_list->sublist("operators") = op_list;

        // insert the flow BC sublist
        Teuchos::ParameterList flow_bc; // = flow_list->sublist("boundary conditions");
        flow_bc = CreateSS_FlowBC_List_(plist);
        flow_list->sublist("boundary conditions") = flow_bc;
	
        // insert sources, if they exist
        Teuchos::ParameterList flow_src; // = flow_list->sublist("source terms");
        flow_src = CreateFlowSrcList_(plist);
        if (flow_src.begin() != flow_src.end()) {
          flow_list->sublist("source terms") = flow_src;
        }

        // do need to add picard list to initialization?
        bool have_picard_params(false);
        Teuchos::ParameterList pic_params;
        if (use_picard_) {
          if (exe_list.isSublist("Numerical Control Parameters")) {
            Teuchos::ParameterList& ncp_list = exe_list.sublist("Numerical Control Parameters");
            if (ncp_list.isSublist("Unstructured Algorithm")) {
              Teuchos::ParameterList& ua_list = ncp_list.sublist("Unstructured Algorithm");
              if (ua_list.isSublist("Steady-State Pseudo-Time Implicit Solver")) {
                have_picard_params = true;
                pic_params = ua_list.sublist("Steady-State Pseudo-Time Implicit Solver");
              }
            }
          }
        }

        // only include a steady state time integrator list if not transient
        Teuchos::ParameterList& ti_mode_list = exe_list.sublist("Time Integration Mode");

        if (! ti_mode_list.isSublist("Transient")) {
          // create sublists for the steady state time integrator
          Teuchos::ParameterList& sti_list = flow_list->sublist("steady state time integrator");

          // error control options
          Teuchos::Array<std::string> err_opts;
          err_opts.push_back(std::string("pressure"));
          sti_list.set<Teuchos::Array<std::string> >("error control options", err_opts);

          // linear solver
          sti_list.set<std::string>("linear solver", ST_SOLVER);
          sti_list.set<std::string>("preconditioner", ST_PRECOND);
          sti_list.set<std::string>("preconditioner enhancement", "none");

          // pressure-lambda constraints
          Teuchos::ParameterList& sti_plamb = sti_list.sublist("pressure-lambda constraints");
          sti_plamb.set<std::string>("method", "projection");
          sti_plamb.set<bool>("inflow krel correction", true);
          sti_plamb.set<std::string>("linear solver", ST_PLAMB_SOLVER);

          // time integration method
          sti_list.set<std::string>("time integration method","BDF1");
          Teuchos::ParameterList& sti_bdf1 = sti_list.sublist("BDF1");

          // use standard timestep controller type
          sti_bdf1.set<std::string>("timestep controller type", ST_TS_CONTROLLER);
          Teuchos::ParameterList& sti_bdf1_std = sti_bdf1.sublist("timestep controller standard parameters");
          sti_bdf1_std.set<int>("max iterations", ST_MAX_ITER);
          sti_bdf1_std.set<int>("min iterations", ST_MIN_ITER);
          sti_bdf1_std.set<double>("time step increase factor", ST_TS_INC_FACTOR);
          sti_bdf1_std.set<double>("time step reduction factor", ST_TS_RED_FACTOR);
          sti_bdf1_std.set<double>("max time step", MAXIMUM_TIME_STEP);
          sti_bdf1_std.set<double>("min time step", ST_MIN_TS);

	  Teuchos::ParameterList* sti_bdf1_solver;

          // solver type
	  if (nonlinear_solver == std::string("Newton") ||
	      nonlinear_solver == std::string("Newton-Picard")) {
	    sti_bdf1.set<std::string>("solver type", "Newton");
	    Teuchos::ParameterList& test = sti_bdf1.sublist("Newton parameters");
	    sti_bdf1_solver = &test;
	    sti_bdf1_solver->set<double>("nonlinear tolerance", STEADY_NONLINEAR_TOLERANCE);
	    sti_bdf1_solver->set<double>("diverged tolerance", ST_NKA_DIVGD_TOL);
	    sti_bdf1_solver->set<double>("max du growth factor", ST_DIVERG_FACT);
	    sti_bdf1_solver->set<int>("max divergent iterations", ST_MAX_DIVERGENT_ITERATIONS);
	    sti_bdf1_solver->set<int>("limit iterations", ST_LIMIT_ITER);
	    sti_bdf1_solver->set<bool>("modify correction", true);
	  }
	  else if (nonlinear_solver == std::string("NKA")) {
	    sti_bdf1.set<std::string>("solver type", "nka");
	    Teuchos::ParameterList& test = sti_bdf1.sublist("nka parameters");
	    sti_bdf1_solver = &test;
	    sti_bdf1_solver->set<double>("nonlinear tolerance", STEADY_NONLINEAR_TOLERANCE);
	    sti_bdf1_solver->set<double>("diverged tolerance", ST_NKA_DIVGD_TOL);
	    sti_bdf1_solver->set<double>("max du growth factor", ST_DIVERG_FACT);
	    sti_bdf1_solver->set<int>("max divergent iterations", ST_MAX_DIVERGENT_ITERATIONS);
	    sti_bdf1_solver->set<int>("max nka vectors", ST_NKA_NUMVEC);
	    sti_bdf1_solver->set<int>("limit iterations", ST_LIMIT_ITER);
	    sti_bdf1_solver->set<bool>("modify correction", modify_correction);
	  }
	  else if (nonlinear_solver == std::string("JFNK")) {
	    sti_bdf1.set<std::string>("solver type", "JFNK");
	    Teuchos::ParameterList& test = sti_bdf1.sublist("JFNK parameters");
	    sti_bdf1_solver = &test;
	    sti_bdf1_solver->set<double>("typical solution value", 1.0);

            Teuchos::ParameterList& tmp = sti_bdf1_solver->sublist("nonlinear solver");     
	    tmp.set<std::string>("solver type", "Newton");
            Teuchos::ParameterList& sti_bdf1_newton = tmp.sublist("Newton parameters");     
	    sti_bdf1_newton.set<double>("diverged tolerance", ST_NKA_DIVGD_TOL);
	    sti_bdf1_newton.set<double>("max du growth factor", ST_DIVERG_FACT);
	    sti_bdf1_newton.set<int>("max divergent iterations", ST_MAX_DIVERGENT_ITERATIONS);
	    sti_bdf1_newton.set<int>("max nka vectors", ST_NKA_NUMVEC);
	    sti_bdf1_newton.set<int>("limit iterations", ST_LIMIT_ITER);

            Teuchos::ParameterList& sti_bdf1_jfmat = sti_bdf1_solver->sublist("JF matrix parameters");
            sti_bdf1_jfmat.set<double>("finite difference epsilon", 1.0e-8);
            sti_bdf1_jfmat.set<std::string>("method for epsilon", "Knoll-Keyes");

            Teuchos::ParameterList& sti_bdf1_linop = sti_bdf1_solver->sublist("linear operator");
            sti_bdf1_linop.set<std::string>("iterative method", "gmres");
            Teuchos::ParameterList& sti_bdf1_gmres = sti_bdf1_linop.sublist("gmres parameters");
            sti_bdf1_gmres.set<double>("error tolerance", 1e-7);
            sti_bdf1_gmres.set<int>("maximum number of iterations", 100);
            std::vector<std::string> criteria;
            criteria.push_back("relative rhs");
            criteria.push_back("relative residual");
            sti_bdf1_gmres.set<Teuchos::Array<std::string> >("convergence criteria", criteria);
          } 
          else {
            msg << "In the definition of Nonlinear Solver Type: you must specify either "
                << "'NKA', 'Newton', 'JFNK', or 'Newton-Picard'";
            Exceptions::amanzi_throw(msg);
          }

          // remaining BDF1 parameters
          sti_bdf1.set<int>("max preconditioner lag iterations", ST_MAX_PREC_LAG);
          sti_bdf1.set<bool>("extrapolate initial guess", true);

          if (exe_list.isSublist("Numerical Control Parameters")) {
            Teuchos::ParameterList& ncp_list = exe_list.sublist("Numerical Control Parameters");
            if (ncp_list.isSublist("Unstructured Algorithm")) {
              Teuchos::ParameterList& ua_list = ncp_list.sublist("Unstructured Algorithm");
              if (ua_list.isSublist("Steady-State Implicit Time Integration")) {
                Teuchos::ParameterList& num_list = ua_list.sublist("Steady-State Implicit Time Integration");
                sti_bdf1_std.set<int>("max iterations", num_list.get<int>("steady max iterations", ST_MAX_ITER));
                sti_bdf1_std.set<int>("min iterations", num_list.get<int>("steady min iterations", ST_MIN_ITER));
                sti_bdf1_solver->set<int>("limit iterations", 
                    num_list.get<int>("steady limit iterations", ST_LIMIT_ITER));
                sti_bdf1_solver->set<double>("nonlinear tolerance",
                    num_list.get<double>("steady nonlinear tolerance", STEADY_NONLINEAR_TOLERANCE));
                sti_bdf1_std.set<double>("time step reduction factor",
                    num_list.get<double>("steady time step reduction factor", ST_TS_RED_FACTOR));
                sti_bdf1_std.set<double>("time step increase factor",
                    num_list.get<double>("steady time step increase factor", ST_TS_INC_FACTOR));
                // sti_bdf1_std.set<double>("max time step",
		//     num_list.get<double>("steady max time step", ST_MAX_TS));
                sti_bdf1.set<int>("max preconditioner lag iterations",
                    num_list.get<int>("steady max preconditioner lag iterations", ST_MAX_PREC_LAG));
                sti_bdf1_solver->set<int>("max divergent iterations",
                    num_list.get<int>("steady max divergent iterations", ST_MAX_DIVERGENT_ITERATIONS));
                sti_bdf1.set<double>("nonlinear iteration damping factor",
                    num_list.get<double>("steady nonlinear iteration damping factor", ST_NONLIN_DAMP));
                sti_bdf1.set<int>("nonlinear iteration initial guess extrapolation order",
                    num_list.get<int>("steady nonlinear iteration initial guess extrapolation order",
                    ST_NONLIN_INIT_GUESS_EXTR_ORD));
                sti_bdf1.set<double>("restart tolerance relaxation factor",
                    num_list.get<double>("steady restart tolerance relaxation factor", ST_NONLIN_INIT_TS_FACTOR));
                sti_bdf1.set<double>("restart tolerance relaxation factor damping",
                    num_list.get<double>("steady restart tolerance relaxation factor damping",
                    ST_NONLIN_INIT_TS_FACTOR_DAMP));
                sti_bdf1_solver->set<double>("max du growth factor",
                    num_list.get<double>("steady nonlinear iteration divergence factor", ST_DIVERG_FACT));

                sti_list.set<Teuchos::Array<std::string> >("error control options",
                    num_list.get<Teuchos::Array<std::string> >("steady error control options",
                    err_opts));
                sti_list.set<std::string>("preconditioner",
                    num_list.get<std::string>("steady preconditioner", ST_PRECOND));

                if (flow_model == "Single Phase") {
                  sti_bdf1_std.set<double>("time step increase factor",
                      num_list.get<double>("steady time step increase factor",ST_SP_DT_INCR_FACTOR));
                }
		// initialization
		if (num_list.get<bool>("steady initialize with darcy", ST_INIT_DARCY_BOOL)) {
		  Teuchos::ParameterList& sti_init = sti_list.sublist("initialization");
                  if (have_picard_params && use_picard_) {
		    sti_init.set<std::string>("method", "picard");
                    sti_init.set<std::string>("linear solver",
                        pic_params.get<std::string>("pseudo time integrator linear solver", PIC_SOLVE));
                    sti_init.set<double>("clipping saturation value", 
                        pic_params.get<double>("pseudo time integrator clipping saturation value", PIC_CLIP_SAT));

                    Teuchos::ParameterList& pic_list = sti_init.sublist("picard parameters");
                    pic_list.set<double>("convergence tolerance",
                        pic_params.get<double>("pseudo time integrator picard convergence tolerance", PICARD_TOLERANCE));
                    pic_list.set<int>("maximum number of iterations",
                        pic_params.get<int>("pseudo time integrator picard maximum number of iterations", PIC_MAX_ITER));
                  } else if (use_picard_) {
		    sti_init.set<std::string>("method", "picard");
                    sti_init.set<double>("clipping saturation value", PIC_CLIP_SAT);

                    Teuchos::ParameterList& pic_list = sti_init.sublist("picard parameters");
                    pic_list.set<std::string>("linear solver", PIC_SOLVE);
                    pic_list.set<double>("convergence tolerance", PICARD_TOLERANCE);
                    pic_list.set<int>("maximum number of iterations", PIC_MAX_ITER);
                  } else {
		    sti_init.set<std::string>("method", "saturated solver");
                    sti_init.set<std::string>("linear solver", ST_INIT_SOLVER);
                  }
		}
              }
            }
          }

          // overwrite parameters for special solvers
          if (nonlinear_solver == std::string("Newton") || 
              nonlinear_solver == std::string("Newton-Picard")) {
            sti_bdf1.set<int>("max preconditioner lag iterations", 0);
	    sti_bdf1.set<bool>("extrapolate initial guess", false);	    
            sti_list.set<std::string>("linear solver", "GMRES for Newton");
            sti_list.set<std::string>("preconditioner enhancement", "GMRES for Newton");
          }

	  if (time_regime == STEADY_REGIME) {
	    flow_list->sublist("time integrator") = sti_list;
	  }
        }


        // only include the transient list if not in steady mode
        if (! ti_mode_list.isSublist("Steady")) {
          // create sublists for the transient state time integrator
          Teuchos::ParameterList& tti_list = flow_list->sublist("transient time integrator");

          // error control options
          Teuchos::Array<std::string> err_opts;
          err_opts.push_back(std::string("pressure"));
          err_opts.push_back(std::string("residual"));
          tti_list.set<Teuchos::Array<std::string> >("error control options", err_opts);

          // linear solver
          if (flow_single_phase) {
            tti_list.set<std::string>("linear solver", TR_SOLVER_DARCY);
          } else if (nonlinear_solver == std::string("Newton")) {
            tti_list.set<std::string>("linear solver", "GMRES for Newton");
          } else {
            tti_list.set<std::string>("linear solver", TR_SOLVER);
          }
          tti_list.set<std::string>("preconditioner", TR_PRECOND);
          tti_list.set<std::string>("preconditioner enhancement", "none");

          // pressure-lambda constraints
          Teuchos::ParameterList& tti_plamb = tti_list.sublist("pressure-lambda constraints");
          tti_plamb.set<std::string>("method", "projection");
          tti_plamb.set<bool>("inflow krel correction", true);
          tti_plamb.set<std::string>("linear solver", TR_PLAMB_SOLVER);

          // time integration method
          tti_list.set<std::string>("time integration method", "BDF1");
          Teuchos::ParameterList& tti_bdf1 = tti_list.sublist("BDF1");

          // use standard timestep controller type
          tti_bdf1.set<std::string>("timestep controller type", TR_TS_CONTROLLER);
          Teuchos::ParameterList& tti_bdf1_std = tti_bdf1.sublist("timestep controller standard parameters");
          tti_bdf1_std.set<int>("max iterations", TR_MAX_ITER);
          tti_bdf1_std.set<int>("min iterations", TR_MIN_ITER);
          tti_bdf1_std.set<double>("time step increase factor", TR_TS_INC_FACTOR);
          tti_bdf1_std.set<double>("time step reduction factor", TR_TS_RED_FACTOR);
          tti_bdf1_std.set<double>("max time step", MAXIMUM_TIME_STEP);
          tti_bdf1_std.set<double>("min time step", TR_MIN_TS);

	  Teuchos::ParameterList* tti_bdf1_solver;

          // solver type
	  if (nonlinear_solver == std::string("Newton") ||
	      nonlinear_solver == std::string("Newton-Picard")) {
	    tti_bdf1.set<std::string>("solver type", "Newton");
	    Teuchos::ParameterList& test = tti_bdf1.sublist("Newton parameters");
	    tti_bdf1_solver = &test;
	    tti_bdf1_solver->set<double>("nonlinear tolerance", TRANSIENT_NONLINEAR_TOLERANCE);
	    tti_bdf1_solver->set<double>("diverged tolerance", TR_NKA_DIVGD_TOL);
	    tti_bdf1_solver->set<double>("max du growth factor", TR_DIVERG_FACT);
	    tti_bdf1_solver->set<int>("max divergent iterations", TR_MAX_DIVERGENT_ITERATIONS);
	    tti_bdf1_solver->set<int>("limit iterations", TR_LIMIT_ITER);
	    tti_bdf1_solver->set<bool>("modify correction", true);
	  }
	  else {
	    tti_bdf1.set<std::string>("solver type", "nka");
	    Teuchos::ParameterList& test = tti_bdf1.sublist("nka parameters");
	    tti_bdf1_solver = &test;
	    tti_bdf1_solver->set<double>("nonlinear tolerance", TRANSIENT_NONLINEAR_TOLERANCE);
	    tti_bdf1_solver->set<double>("diverged tolerance", TR_NKA_DIVGD_TOL);
	    tti_bdf1_solver->set<double>("max du growth factor", TR_DIVERG_FACT);
	    tti_bdf1_solver->set<int>("max divergent iterations", TR_MAX_DIVERGENT_ITERATIONS);
	    tti_bdf1_solver->set<int>("max nka vectors", TR_NKA_NUMVEC);
	    tti_bdf1_solver->set<int>("limit iterations", TR_LIMIT_ITER);
	    tti_bdf1_solver->set<bool>("modify correction", modify_correction);
	  }

          // remaining parameters
          tti_bdf1.set<int>("max preconditioner lag iterations", TR_MAX_PREC_LAG);
          tti_bdf1.set<bool>("extrapolate initial guess", true);

          if (exe_list.isSublist("Numerical Control Parameters")) {
            Teuchos::ParameterList& ncp_list = exe_list.sublist("Numerical Control Parameters");
            if (ncp_list.isSublist("Unstructured Algorithm")) {
              Teuchos::ParameterList& ua_list = ncp_list.sublist("Unstructured Algorithm");
              if (ua_list.isSublist("Transient Implicit Time Integration")) {
                Teuchos::ParameterList& num_list = ua_list.sublist("Transient Implicit Time Integration");

                tti_bdf1_std.set<int>("max iterations",
                    num_list.get<int>("transient max iterations", TR_MAX_ITER));
                tti_bdf1_std.set<int>("min iterations",
                    num_list.get<int>("transient min iterations", TR_MIN_ITER));
                tti_bdf1_solver->set<int>("limit iterations", 
                    num_list.get<int>("transient limit iterations", TR_LIMIT_ITER));
                tti_bdf1_solver->set<double>("nonlinear tolerance",
                    num_list.get<double>("transient nonlinear tolerance", TRANSIENT_NONLINEAR_TOLERANCE));
                tti_bdf1_std.set<double>("time step reduction factor",
                    num_list.get<double>("transient time step reduction factor", TR_TS_RED_FACTOR));
                tti_bdf1_std.set<double>("time step increase factor",
                    num_list.get<double>("transient time step increase factor", TR_TS_INC_FACTOR));
                // tti_bdf1_std.set<double>("max time step", num_list.get<double>("transient max time step", TR_MAX_TS));
                tti_bdf1.set<int>("max preconditioner lag iterations",
                    num_list.get<int>("transient max preconditioner lag iterations", TR_MAX_PREC_LAG));
                tti_bdf1_solver->set<int>("max divergent iterations",
                    num_list.get<int>("transient max divergent iterations", TR_MAX_DIVERGENT_ITERATIONS));
                tti_bdf1.set<double>("nonlinear iteration damping factor",
                    num_list.get<double>("transient nonlinear iteration damping factor", TR_NONLIN_DAMP));
                tti_bdf1.set<int>("nonlinear iteration initial guess extrapolation order",
                    num_list.get<int>("transient nonlinear iteration initial guess extrapolation order", 
                    TR_NONLIN_INIT_GUESS_EXTR_ORD));
                tti_bdf1.set<double>("restart tolerance relaxation factor",
                    num_list.get<double>("transient restart tolerance relaxation factor", TR_NONLIN_INIT_TS_FACTOR));
                tti_bdf1.set<double>("restart tolerance relaxation factor damping",
                    num_list.get<double>("transient restart tolerance relaxation factor damping", 
                    TR_NONLIN_INIT_TS_FACTOR_DAMP));
                tti_bdf1_solver->set<double>("max du growth factor",
                    num_list.get<double>("transient nonlinear iteration divergence factor", TR_DIVERG_FACT));

                tti_list.set<Teuchos::Array<std::string> >("error control options",
                    num_list.get<Teuchos::Array<std::string> >("transient error control options",
                    err_opts));
                tti_list.set<std::string>("preconditioner",
                    num_list.get<std::string>("transient preconditioner", TR_PRECOND));

                if (flow_model == "Single Phase") {
                  tti_bdf1_std.set<double>("time step increase factor",
                      num_list.get<double>("transient time step increase factor", TR_SP_DT_INCR_FACTOR));
                }
		// create an initialization sublist
		if (num_list.get<bool>("transient initialize with darcy", TR_INIT_DARCY_BOOL)) {
		  Teuchos::ParameterList& tti_init = tti_list.sublist("initialization");
		  tti_init.set<std::string>("method", "saturated solver");
		  tti_init.set<std::string>("linear solver", TR_INIT_SOLVER);
		}
	      }
            }
          }

          if (nonlinear_solver == std::string("Newton")) {
            tti_bdf1.set<int>("max preconditioner lag iterations", 0);
	    tti_bdf1.set<bool>("extrapolate initial guess", false);
            tti_list.set<std::string>("preconditioner enhancement", "GMRES for Newton");
          }

          tti_list.sublist("VerboseObject") = CreateVerbosityList_(verbosity_level);

	  if (time_regime == TRANSIENT_REGIME) { 
	    flow_list->sublist("time integrator") = tti_list; 
	  }
	}
      }
    }
  }

  // cleaning flow for better testing of PKs
  flow_list->remove("steady state time integrator", false);
  flow_list->remove("transient time integrator", false);

  return flw_list;
 }


/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateFlowSrcList_(Teuchos::ParameterList* plist)
{
  Errors::Message msg;
  Teuchos::ParameterList src_list;
  Teuchos::ParameterList& src_sublist = plist->sublist("Sources");

  for (Teuchos::ParameterList::ConstIterator i = src_sublist.begin(); i != src_sublist.end(); i++) {
    // look at sublists
    if (src_sublist.isSublist(src_sublist.name(i))) {
      Teuchos::ParameterList& src = src_sublist.sublist(src_sublist.name(i));
      Teuchos::ParameterList src_fn;
      Teuchos::ParameterList src_sub_out;      

      if (src.isSublist("Source: Volume Weighted")) {
        src_sub_out.set<std::string>("spatial distribution method","volume");
        src_fn = src.sublist("Source: Volume Weighted");
	if (!src.sublist("Source: Volume Weighted").isParameter("Values")) continue;
      }
      else if (src.isSublist("Source: Permeability Weighted")) {
        src_sub_out.set<std::string>("spatial distribution method","permeability");
        src_fn = src.sublist("Source: Permeability Weighted");
	if (!src.sublist("Source: Permeability Weighted").isParameter("Values")) continue;
      }
      else if (src.isSublist("Source: Uniform")) {
        src_sub_out.set<std::string>("spatial distribution method","none");
        src_fn = src.sublist("Source: Uniform");
	if (!src.sublist("Source: Uniform").isParameter("Values")) continue;
      }
      else {
        msg << "In the definition of Sources: you must either specify "
            << "'Source: Volume Weighted', 'Source: Permeability Weighted', or 'Source: Uniform'";
        Exceptions::amanzi_throw(msg);
      }

      std::string name = src_sublist.name(i);

      // get the regions
      Teuchos::Array<std::string> regions = src.get<Teuchos::Array<std::string> >("Assigned Regions");
      src_sub_out.set<Teuchos::Array<std::string> >("regions", regions);
      vv_src_regions.insert(vv_src_regions.end(), regions.begin(), regions.end());

      // create time function
      Teuchos::ParameterList& src_sub_out_fn = src_sub_out.sublist("sink");
      Teuchos::Array<double> values = src_fn.get<Teuchos::Array<double> >("Values");

      // write the native time function
      if (values.size() == 1) {
        src_sub_out_fn.sublist("function-constant").set<double>("value", values[0]);
      } 
      else if (values.size() > 1) {
        Teuchos::Array<double> times = src_fn.get<Teuchos::Array<double> >("Times");
        Teuchos::Array<std::string> time_fns = src_fn.get<Teuchos::Array<std::string> >("Time Functions");

        Teuchos::ParameterList& ssofn = src_sub_out_fn.sublist("function-tabular");

        ssofn.set<Teuchos::Array<double> >("x values", times);
        ssofn.set<Teuchos::Array<double> >("y values", values);

        Teuchos::Array<std::string> forms_(time_fns.size());

        for (int i = 0; i < time_fns.size(); i++) {
          if (time_fns[i] == "Linear") {
            forms_[i] = "linear";
          }
          else if (time_fns[i] == "Constant") {
            forms_[i] = "constant";
          }
          else {
            msg << "In the definition of Sources: time function can only be 'Linear' or 'Constant'";
            Exceptions::amanzi_throw(msg);
          }
        }
        ssofn.set<Teuchos::Array<std::string> >("forms",forms_);
      }
      else {
        msg << "In the definition of Sources: something is wrong with the input";
        Exceptions::amanzi_throw(msg);
      }
      src_list.sublist(name) = src_sub_out;
    }
  }

  return src_list;
}


/* ******************************************************************
* DPC sublist
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateSS_FlowBC_List_(Teuchos::ParameterList* plist)
{
  Errors::Message msg;
  Teuchos::ParameterList ssf_list;
  Teuchos::ParameterList& bc_sublist = plist->sublist("Boundary Conditions");

  int bc_counter = 0;

  for (Teuchos::ParameterList::ConstIterator i = bc_sublist.begin(); i != bc_sublist.end(); i++) {
    // look at sublists
    if (bc_sublist.isSublist(bc_sublist.name(i))) {
      Teuchos::ParameterList& bc = bc_sublist.sublist(bc_sublist.name(i));

      // get the regions
      Teuchos::Array<std::string> regions = bc.get<Teuchos::Array<std::string> >("Assigned Regions");
      vv_bc_regions.insert(vv_bc_regions.end(), regions.begin(), regions.end());

      if (bc.isSublist("BC:Zero Flow")) {
        // this is the natural BC for flow and we need not list it explicitly
      }

      else if (bc.isSublist("BC: Flux")) {
        Teuchos::ParameterList& bc_flux = bc.sublist("BC: Flux");

        Teuchos::Array<double> times = bc_flux.get<Teuchos::Array<double> >("Times");
        Teuchos::Array<std::string> time_fns = bc_flux.get<Teuchos::Array<std::string> >("Time Functions");

        Teuchos::Array<double> flux;

        if (bc_flux.isParameter("Inward Mass Flux")) {
          flux = bc_flux.get<Teuchos::Array<double> >("Inward Mass Flux");
          for (int i = 0; i < flux.size(); i++) flux[i] *= -1;
        } else if (bc_flux.isParameter("Outward Mass Flux")) {
          flux = bc_flux.get<Teuchos::Array<double> >("Outward Mass Flux");
        } else if (bc_flux.isParameter("Inward Volumetric Flux")) {
          flux = bc_flux.get<Teuchos::Array<double> >("Inward Volumetric Flux");
          for (int i = 0; i < flux.size(); i++) flux[i] *= -constant_density;
        } else if (bc_flux.isParameter("Outward Volumetric Flux")) {
          flux = bc_flux.get<Teuchos::Array<double> >("Outward Volumetric Flux");
          for (int i = 0; i < flux.size(); i++) flux[i] *= constant_density;
        }

        std::stringstream ss;
        ss << "BC " << bc_counter++;

        Teuchos::ParameterList& tbc = ssf_list.sublist("mass flux").sublist(ss.str());
        tbc.set<Teuchos::Array<std::string> >("regions", regions );

        tbc.set<bool>("rainfall", bc_flux.get<bool>("rainfall",false));

        if (times.size() == 1) {
          Teuchos::ParameterList& tbcs = tbc.sublist("outward mass flux").sublist("function-constant");
          tbcs.set<double>("value",flux[0]);
        } else {
          Teuchos::ParameterList& tbcs = tbc.sublist("outward mass flux").sublist("function-tabular");

          tbcs.set<Teuchos::Array<double> >("x values", times);
          tbcs.set<Teuchos::Array<double> >("y values", flux);

          std::vector<std::string> forms_(time_fns.size());

          for (int i = 0; i < time_fns.size(); i++) {

            if (time_fns[i] == "Linear") {
              forms_[i] = "linear";
            } else if (time_fns[i] == "Constant") {
              forms_[i] = "constant";
            } else {
              msg << "In the definition of BCs: tabular function can only be Linear or Constant";
              Exceptions::amanzi_throw(msg);
            }
          }

          Teuchos::Array<std::string> forms = forms_;
          tbcs.set<Teuchos::Array<std::string> >("forms", forms);
        }
      }
      else if (bc.isSublist("BC: Uniform Pressure")) {
        Teuchos::ParameterList& bc_dir = bc.sublist("BC: Uniform Pressure");

        Teuchos::Array<double>      times = bc_dir.get<Teuchos::Array<double> >("Times");
        Teuchos::Array<std::string> time_fns = bc_dir.get<Teuchos::Array<std::string> >("Time Functions");
        Teuchos::Array<double>      values = bc_dir.get<Teuchos::Array<double> >("Values");

        std::stringstream ss;
        ss << "BC " << bc_counter++;

        Teuchos::ParameterList& tbc = ssf_list.sublist("pressure").sublist(ss.str());
        tbc.set<Teuchos::Array<std::string> >("regions", regions);

        if (times.size() == 1) {
          Teuchos::ParameterList& tbcs = tbc.sublist("boundary pressure").sublist("function-constant");
          tbcs.set<double>("value", values[0]);
        } else {
          Teuchos::ParameterList& tbcs = tbc.sublist("boundary pressure").sublist("function-tabular");

          tbcs.set<Teuchos::Array<double> >("x values", times);
          tbcs.set<Teuchos::Array<double> >("y values", values);

          std::vector<std::string> forms_(time_fns.size());

          for (int i = 0; i < time_fns.size(); i++)
            if (time_fns[i] == "Linear") {
              forms_[i] = "linear";
            } else if (time_fns[i] == "Constant") {
              forms_[i] = "constant";
            } else {
              Exceptions::amanzi_throw(Errors::Message("Tabular function can only be Linear or Constant"));
            }
          Teuchos::Array<std::string> forms = forms_;
          tbcs.set<Teuchos::Array<std::string> >("forms", forms);
        }

      }
      else if (bc.isSublist("BC: Linear Pressure")) {
        Teuchos::ParameterList& bc_dir = bc.sublist("BC: Linear Pressure");
        Teuchos::Array<double> grad = bc_dir.get<Teuchos::Array<double> >("Gradient Value");
        Teuchos::Array<double> refcoord = bc_dir.get<Teuchos::Array<double> >("Reference Point");
        double refval = bc_dir.get<double>("Reference Value");

        Teuchos::Array<double> grad_with_time(grad.size() + 1);
        grad_with_time[0] = 0.0;
        for (int j = 0; j != grad.size(); ++j) {
          grad_with_time[j + 1] = grad[j];
        }

        Teuchos::Array<double> refcoord_with_time(refcoord.size() + 1);
        refcoord_with_time[0] = 0.0;
        for (int j = 0; j != refcoord.size(); ++j) {
          refcoord_with_time[j + 1] = refcoord[j];
        }

        std::stringstream ss;
        ss << "BC " << bc_counter++;

        Teuchos::ParameterList& tbc = ssf_list.sublist("pressure").sublist(ss.str());
        tbc.set<Teuchos::Array<std::string> >("regions", regions);

        Teuchos::ParameterList& tbcs = tbc.sublist("boundary pressure").sublist("function-linear");
        tbcs.set<double>("y0", refval);
        tbcs.set<Teuchos::Array<double> >("x0", refcoord_with_time);
        tbcs.set<Teuchos::Array<double> >("gradient", grad_with_time);

      }
      else if (bc.isSublist("BC: Linear Hydrostatic")) {
        Teuchos::ParameterList& bc_dir = bc.sublist("BC: Linear Hydrostatic");
        Teuchos::Array<double> grad = bc_dir.get<Teuchos::Array<double> >("Gradient Value");
        Teuchos::Array<double> refcoord = bc_dir.get<Teuchos::Array<double> >("Reference Point");
        double refval = bc_dir.get<double>("Reference Water Table Height");
        std::string submodel = bc_dir.get<std::string>("Submodel", "None");

        Teuchos::Array<double> grad_with_time(grad.size() + 1);
        grad_with_time[0] = 0.0;
        for (int j = 0; j != grad.size(); ++j) {
          grad_with_time[j + 1] = grad[j];
        }

        Teuchos::Array<double> refcoord_with_time(refcoord.size() + 1);
        refcoord_with_time[0] = 0.0;
        for (int j = 0; j != refcoord.size(); ++j) {
          refcoord_with_time[j + 1] = refcoord[j];
        }

        std::stringstream ss;
        ss << "BC " << bc_counter++;

        Teuchos::ParameterList& tbc = ssf_list.sublist("static head").sublist(ss.str());
        tbc.set<Teuchos::Array<std::string> >("regions", regions);

        Teuchos::ParameterList& tbcs = tbc.sublist("water table elevation").sublist("function-linear");
        tbcs.set<double>("y0", refval);
        tbcs.set<Teuchos::Array<double> >("x0", refcoord_with_time);
        tbcs.set<Teuchos::Array<double> >("gradient", grad_with_time);

        if (submodel == "No Flow Above Water Table") {
          tbc.set<bool>("no flow above water table", true);
        } else if (submodel != "None") {
          msg << "In 'BC: Linear Hydrostatic': optional parameter 'Submodel': valid values are 'No Flow"
              << " Above Water Table' or 'None'";
          Exceptions::amanzi_throw(msg);
        } 
      }
      else if (bc.isSublist("BC: Hydrostatic")) {
        Teuchos::ParameterList& bc_dir = bc.sublist("BC: Hydrostatic");

        Teuchos::Array<double> times = bc_dir.get<Teuchos::Array<double> >("Times");
        Teuchos::Array<std::string> time_fns = bc_dir.get<Teuchos::Array<std::string> >("Time Functions");
        Teuchos::Array<double> values = bc_dir.get<Teuchos::Array<double> >("Water Table Height");
        std::string coordsys = bc_dir.get<std::string>("Coordinate System", BCHYDRST_COORD);
        std::string submodel = bc_dir.get<std::string>("Submodel", "None");

        std::stringstream ss;
        ss << "BC " << bc_counter++;

        Teuchos::ParameterList& tbc = ssf_list.sublist("static head").sublist(ss.str());
        tbc.set<Teuchos::Array<std::string> >("regions", regions );

        if (coordsys == "Absolute") {
          tbc.set<bool>("relative to top", false);
        } else if (coordsys == "Relative") {
          tbc.set<bool>("relative to top", true);
        } else {
          // we have a default for this value... "Absolute", if for some reason this does not
          // get read, then we must bail
          msg << "In 'BC: Hydrostatic': must specify a value for parameter 'Coordinate System',"
              << " valid values are 'Absolute' and 'Relative'";
          Exceptions::amanzi_throw(msg);
        }

        if (times.size() == 1) {
          Teuchos::ParameterList& tbcs = tbc.sublist("water table elevation").sublist("function-constant");
          tbcs.set<double>("value",values[0]);
        } else {
          Teuchos::ParameterList& tbcs = tbc.sublist("water table elevation").sublist("function-tabular");

          tbcs.set<Teuchos::Array<double> >("x values", times);
          tbcs.set<Teuchos::Array<double> >("y values", values);

          std::vector<std::string> forms_(time_fns.size());

          for (int i = 0; i < time_fns.size(); i++)
            if (time_fns[i] == "Linear") {
              forms_[i] = "linear";
            } else if (time_fns[i] == "Constant") {
              forms_[i] = "constant";
            } else  {
              Exceptions::amanzi_throw(Errors::Message("Tabular function can only be Linear or Constant"));
            }
          Teuchos::Array<std::string> forms = forms_;
          tbcs.set<Teuchos::Array<std::string> >("forms", forms);
        }

        if (submodel == "No Flow Above Water Table") {
          tbc.set<bool>("no flow above water table", true);
        } else if (submodel != "None") {
          msg << "In 'BC: Hydrostatic': optional parameter 'Submodel': valid values are 'No Flow Above"
              << " Water Table' or 'None'";
          Exceptions::amanzi_throw(msg);
        } 
      }
      else if (bc.isSublist("BC: Seepage")) {
        Teuchos::ParameterList& bc_flux = bc.sublist("BC: Seepage");

        Teuchos::Array<double> times = bc_flux.get<Teuchos::Array<double> >("Times");
        Teuchos::Array<std::string> time_fns = bc_flux.get<Teuchos::Array<std::string> >("Time Functions");

        Teuchos::Array<double> flux;

        if (bc_flux.isParameter("Inward Mass Flux")) {
          flux = bc_flux.get<Teuchos::Array<double> >("Inward Mass Flux");
          for (int i = 0; i < flux.size(); i++) flux[i] *= -1;
        } else if (bc_flux.isParameter("Inward Volumetric Flux")) {
          for (int i = 0; i < flux.size(); i++) flux[i] *= -constant_density;
        } else {
          msg << "In \"BC: Seepage\" we can only handle \"Inward Mass Flux\" or \"Inward Volumetric Flux\"";
          Exceptions::amanzi_throw(msg);
        }

        std::stringstream ss;
        ss << "BC " << bc_counter++;

        Teuchos::ParameterList& tbc = ssf_list.sublist("seepage face").sublist(ss.str());
        tbc.set<Teuchos::Array<std::string> >("regions", regions );
        tbc.set<bool>("rainfall", bc_flux.get<bool>("rainfall",false));

        if (times.size() == 1) {
          Teuchos::ParameterList& tbcs = tbc.sublist("outward mass flux").sublist("function-constant");
          tbcs.set<double>("value",flux[0]);
        } else {
          Teuchos::ParameterList& tbcs = tbc.sublist("outward mass flux").sublist("function-tabular");

          tbcs.set<Teuchos::Array<double> >("x values", times);
          tbcs.set<Teuchos::Array<double> >("y values", flux);

          std::vector<std::string> forms_(time_fns.size());

          for (int i = 0; i < time_fns.size(); i++) {
            if (time_fns[i] == "Linear") {
              forms_[i] = "linear";
            } else if (time_fns[i] == "Constant") {
              forms_[i] = "constant";
            } else {
              msg << "In the definition of BCs: tabular function can only be Linear or Constant";
              Exceptions::amanzi_throw(msg);
            }
          }

          Teuchos::Array<std::string> forms = forms_;
          tbcs.set<Teuchos::Array<std::string> >("forms", forms);
        }

        // hack: select one region for transport diagnostics
        transport_diagnostics_.insert(transport_diagnostics_.end(), regions.begin(), regions.end());
      }
      // TODO...
      // add the rest of the boundary conditions
    }
  }

  return ssf_list;
}


/* ******************************************************************
* WRM sublist
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateWRM_List_(Teuchos::ParameterList* plist)
{
  Errors::Message msg;
  Teuchos::OSTab tab = vo_->getOSTab();

  Teuchos::ParameterList wrm_list;

  // loop through the material properties list and extract the water retention model info
  Teuchos::ParameterList& matprop_list = plist->sublist("Material Properties");

  int counter = 0;
  for (Teuchos::ParameterList::ConstIterator i = matprop_list.begin(); i != matprop_list.end(); i++) {
    // get the wrm parameters
    Teuchos::ParameterList& cp_list = matprop_list.sublist(i->first);

    // we can have either van Genuchten or Brooks Corey
    if (cp_list.isSublist("Capillary Pressure: van Genuchten")) {
      Teuchos::ParameterList vG_list = cp_list.sublist("Capillary Pressure: van Genuchten");
      std::string rel_perm = vG_list.get<std::string>("Relative Permeability");
      if (rel_perm != "Mualem" && rel_perm != "Burdine") {
        std::stringstream ss;
        ss << "Currently we only have Mualem or Burdine as the relative permeability models";
        Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
      }

      double alpha = vG_list.get<double>("alpha");
      double Sr = vG_list.get<double>("Sr");
      double m = vG_list.get<double>("m");
      double ell;
      if (rel_perm == "Mualem") {
        ell = vG_list.get<double>("ell",ELL_MUALEM);
      } else if (rel_perm == "Burdine") {
        ell = vG_list.get<double>("ell",ELL_BURDINE);
      }
      double krel_smooth = vG_list.get<double>("krel smoothing interval", 0.0);
      if (krel_smooth < 0.0) {
        Exceptions::amanzi_throw(Errors::Message("If krel smoothing interval is specified it must be positive."));
      }

      // now get the assigned regions
      Teuchos::Array<std::string> regions = cp_list.get<Teuchos::Array<std::string> >("Assigned Regions");

      for (Teuchos::Array<std::string>::const_iterator i = regions.begin(); i != regions.end(); i++) {
        std::stringstream ss;
        ss << "WRM for " << *i;

        Teuchos::ParameterList& wrm_sublist = wrm_list.sublist(ss.str());

        wrm_sublist.set<std::string>("water retention model", "van Genuchten");
        wrm_sublist.set<std::string>("region", *i);
        wrm_sublist.set<double>("van Genuchten m", m);
        wrm_sublist.set<double>("van Genuchten l", ell);
        wrm_sublist.set<double>("van Genuchten alpha", alpha);
        wrm_sublist.set<double>("residual saturation", Sr);
        wrm_sublist.set<double>("regularization interval", krel_smooth);
        wrm_sublist.set<std::string>("relative permeability model", rel_perm);
      
        if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
          Teuchos::ParameterList& file_list = wrm_sublist.sublist("output");
          std::stringstream name;
          name << *i << ".txt";
          file_list.set<std::string>("file", name.str());
          file_list.set<int>("number of points", 1000);

          *vo_->os() << "water retention curve file:" << name.str() << std::endl;
        }
      }
    } else if (cp_list.isSublist("Capillary Pressure: Brooks Corey")) {
      Teuchos::ParameterList& BC_list = cp_list.sublist("Capillary Pressure: Brooks Corey");

      std::string rel_perm = BC_list.get<std::string>("Relative Permeability");
      if (rel_perm != "Mualem" && rel_perm != "Burdine") {
        std::stringstream ss;
        ss << "Currently we only have Mualem or Burdine as the relative permeability models";
        Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
      }

      double lambda = BC_list.get<double>("lambda");
      double alpha = BC_list.get<double>("alpha");
      double ell;
      if (rel_perm == "Mualem") {
        ell = BC_list.get<double>("ell",ELL_MUALEM);
      } else if (rel_perm == "Burdine") {
        ell = BC_list.get<double>("ell",ELL_BURDINE);
      }
      double Sr = BC_list.get<double>("Sr");
      double krel_smooth = BC_list.get<double>("krel smoothing interval",0.0);

      if (krel_smooth < 0.0) {
        Exceptions::amanzi_throw(Errors::Message("If krel smoothing interval is specified it must be positive."));
      }

      // now get the assigned regions
      Teuchos::Array<std::string> regions = cp_list.get<Teuchos::Array<std::string> >("Assigned Regions");

      for (Teuchos::Array<std::string>::const_iterator i = regions.begin();
           i != regions.end(); i++) {
        std::stringstream ss;
        ss << "WRM for " << *i;

        Teuchos::ParameterList& wrm_sublist = wrm_list.sublist(ss.str());

        wrm_sublist.set<std::string>("water retention model", "Brooks Corey");
        wrm_sublist.set<std::string>("region", *i);
        wrm_sublist.set<double>("Brooks Corey lambda", lambda);
        wrm_sublist.set<double>("Brooks Corey alpha", alpha);
        wrm_sublist.set<double>("Brooks Corey l", ell);
        wrm_sublist.set<double>("residual saturation", Sr);
        wrm_sublist.set<double>("regularization interval", krel_smooth);
        wrm_sublist.set<std::string>("relative permeability model", rel_perm);

        if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
          Teuchos::ParameterList& file_list = wrm_sublist.sublist("output");
          std::stringstream name;
          name << *i << ".txt";
          file_list.set<std::string>("file", name.str());
          file_list.set<int>("number of points", 1000);

          *vo_->os() << "water retention curve file:" << name.str() << std::endl;
        }
      }
    } else {
      msg << "An unknown capillary pressure model was specified, must specify" 
          << " either van Genuchten or Brooks Corey";
      Exceptions::amanzi_throw(msg);
    }
  }
  return wrm_list;
}


/* ******************************************************************
* Flow operators sublist
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateFlowOperatorList_(
    const std::string& disc_method, const std::string& prec_method,
    const std::string& nonlinear_solver, const std::string& rel_perm)
{
  Teuchos::ParameterList op_list;

  Teuchos::ParameterList tmp_list;
  tmp_list.set<std::string>("discretization primary", boost::to_lower_copy(disc_method));
  tmp_list.set<std::string>("discretization secondary", "mfd: optimized for sparsity");

  if (disc_method != "FV: Default"){
    Teuchos::Array<std::string> stensil(2);
    stensil[0] = "face";
    stensil[1] = "cell";
    tmp_list.set<Teuchos::Array<std::string> >("schema", stensil);

    if (prec_method != "Linearized Operator") stensil.remove(1);
    tmp_list.set<Teuchos::Array<std::string> >("preconditioner schema", stensil);
    tmp_list.set<bool>("gravity", true);
  }
  else {
    Teuchos::Array<std::string> stensil(1);
    stensil[0] = "cell";
    tmp_list.set<Teuchos::Array<std::string> >("schema", stensil);

    tmp_list.set<Teuchos::Array<std::string> >("preconditioner schema", stensil);
    tmp_list.set<bool>("gravity", true);
  }

  // create two operators for matrix and preconditioner
  op_list.sublist("diffusion operator").sublist("matrix") = tmp_list;
  op_list.sublist("diffusion operator").sublist("preconditioner") = tmp_list;

  if (prec_method == "Linearized Operator") {
    op_list.sublist("diffusion operator").sublist("preconditioner")
        .set<std::string>("newton correction", "approximate jacobian");
  }

  if (nonlinear_solver == "Newton") {
    Teuchos::ParameterList& prec_list = 
        op_list.sublist("diffusion operator").sublist("preconditioner");
    prec_list.set<std::string>("newton correction", "true jacobian");

    Teuchos::ParameterList& slist = prec_list.sublist("linear operator");
    slist.set<std::string>("iterative method", "gmres");
    Teuchos::ParameterList& gmres_list = slist.sublist("gmres parameters");
    gmres_list.set<double>("error tolerance", TRANSIENT_NONLINEAR_TOLERANCE * 1e-2);
    gmres_list.set<int>("maximum number of iterations", 50);
    std::vector<std::string> criteria;
    criteria.push_back("relative rhs");
    criteria.push_back("relative residual");
    gmres_list.set<Teuchos::Array<std::string> >("convergence criteria", criteria);
    gmres_list.sublist("VerboseObject") = CreateVerbosityList_("low");
  }

  return op_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
