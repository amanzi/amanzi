#include <sstream>
#include <string>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputParserIS.hh"
#include "InputParserIS_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

/* ******************************************************************
* Collects default preconditioners
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreatePreconditionersList_(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList prec_list;
  prec_list.sublist("Trilinos ML") = CreateDPC_List_(plist);
  prec_list.sublist("Hypre AMG") = CreateHypreAMG_List_(plist);
  prec_list.sublist("Block ILU") = CreateBILU_List_(plist);
  return prec_list;
}


/* ******************************************************************
* Collects linear solvers
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateSolversList_(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList solver_list;
  Teuchos::ParameterList& aztecoo_list = solver_list.sublist("AztecOO");
  Teuchos::ParameterList& pcg_list = solver_list.sublist("PCG with Hypre AMG");
  Teuchos::ParameterList& gmres_list = solver_list.sublist("GMRES for Newton");

  // define defaults...
  double tol = LIN_SOLVE_TOL;
  int maxiter = LIN_SOLVE_MAXITER;
  std::string method = LIN_SOLVE_METHOD;
  std::string prec = LIN_SOLVE_PREC;

  // get values from Execution control list if they exist
  if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
    Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters");
    if (ncp_list.isSublist("Unstructured Algorithm")) {
      Teuchos::ParameterList& ncpu_list = ncp_list.sublist("Unstructured Algorithm");
      if (ncpu_list.isSublist("Linear Solver")) {
        Teuchos::ParameterList& num_list = ncpu_list.sublist("Linear Solver");
        if (num_list.isParameter("linear solver tolerance"))
          tol = num_list.get<double>("linear solver tolerance");
        if (num_list.isParameter("linear solver maximum iterations"))
          maxiter = num_list.get<int>("linear solver maximum iterations");
        if (num_list.isParameter("linear solver iterative method"))
          method = num_list.get<std::string>("linear solver iterative method");
        if (num_list.isParameter("linear solver method"))
          prec = num_list.get<std::string>("linear solver preconditioner");
      }
    }
  }
  aztecoo_list.set<std::string>("preconditioner", prec);
  aztecoo_list.set<std::string>("iterative method", method);
  {
    method.append(" parameters");
    Teuchos::ParameterList& method_list = aztecoo_list.sublist(method);
    method_list.set<double>("error tolerance", tol);
    method_list.set<int>("maximum number of iterations", maxiter);
    method_list.sublist("VerboseObject") = CreateVerbosityList_(verbosity_level);
  }

  // add default PCG solver
  pcg_list.set<std::string>("preconditioner", prec);
  pcg_list.set<std::string>("iterative method", "pcg");
  {
    Teuchos::ParameterList& method_list = pcg_list.sublist("pcg parameters");
    method_list.set<double>("error tolerance", tol);
    method_list.set<int>("maximum number of iterations", maxiter);
    method_list.sublist("VerboseObject") = CreateVerbosityList_(verbosity_level);
  }

  // add default "GMRES for Newton" solver
  gmres_list.set<std::string>("iterative method", "gmres");
  {
    Teuchos::ParameterList& method_list = gmres_list.sublist("gmres parameters");
    method_list.set<double>("error tolerance", 1e-7);
    method_list.set<int>("maximum number of iterations", 50);
    std::vector<std::string> criteria;
    criteria.push_back("relative rhs");
    criteria.push_back("relative residual");
    method_list.set<Teuchos::Array<std::string> >("convergence criteria", criteria);
    method_list.sublist("VerboseObject") = CreateVerbosityList_("low");
  }

  return solver_list;
}


/* ******************************************************************
* ML preconditioner sublist
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateDPC_List_(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList dpc_list;

  dpc_list.set<std::string>("preconditioner type", "ml");

  double aggthr(ML_AGG_THR);
  std::string smthtyp(ML_SMOOTHER);
  int ncycles(ML_NCYC);
  int nsmooth(ML_NSMOOTH);
  std::string nonlinear_solver("NKA");

  if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
    Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters");
    if (ncp_list.isSublist("Unstructured Algorithm")) {
      Teuchos::ParameterList& ncpu_list = ncp_list.sublist("Unstructured Algorithm");
      if (ncpu_list.isSublist("Preconditioners")) {
        Teuchos::ParameterList& ncpup_list = ncpu_list.sublist("Preconditioners");
        if (ncpup_list.isSublist("Trilinos ML")) {
          Teuchos::ParameterList& ml_list = ncpup_list.sublist("Trilinos ML");
          if (ml_list.isParameter("ML aggregation threshold")) {
            aggthr = ml_list.get<double>("ML aggregation threshold");
          }
          if (ml_list.isParameter("ML smoother type")) {
            smthtyp = ml_list.get<std::string>("ML smoother type");
          }
          if (ml_list.isParameter("ML cycle applications")) {
            ncycles = ml_list.get<int>("ML cycle applications");
          }
          if (ml_list.isParameter("ML smoother sweeps")) {
            nsmooth = ml_list.get<int>("ML smoother sweeps");
          }
        }
      }
      Teuchos::ParameterList& solver_list = ncpu_list.sublist("Nonlinear Solver");
      if (solver_list.isParameter("Nonlinear Solver Type")) {
	nonlinear_solver = solver_list.get<std::string>("Nonlinear Solver Type");
      }
      if (nonlinear_solver == std::string("Newton")) {
	dpc_list.set<std::string>("discretization method", "fv: default");
      }
      else{
	dpc_list.set<std::string>("discretization method", "mfd: optimized for sparsity");
      }
    }
  }

  Teuchos::ParameterList& ml_list = dpc_list.sublist("ml parameters");
  ml_list.set<int>("ML output", ML_OUTPUT);
  ml_list.set<int>("max levels", ML_MAXLVLS);
  ml_list.set<std::string>("prec type",ML_PRECTYPE);
  ml_list.set<int>("cycle applications", ncycles);
  ml_list.set<std::string>("aggregation: type", ML_AGGTYPE);
  ml_list.set<double>("aggregation: damping factor", ML_AGGDAMP);
  ml_list.set<double>("aggregation: threshold", aggthr);
  ml_list.set<std::string>("eigen-analysis: type",ML_EIGENANAL_TYPE);
  ml_list.set<int>("eigen-analysis: iterations", ML_EIGENANAL_ITERS);
  ml_list.set<int>("smoother: sweeps", nsmooth);
  ml_list.set<double>("smoother: damping factor", ML_SMOOTH_DAMP);
  ml_list.set<std::string>("smoother: pre or post", ML_SMOOTH_PRE_POST);
  ml_list.set<std::string>("smoother: type", smthtyp);
  ml_list.set<double>("smoother: damping factor", ML_SMOOTH_DAMP);
  ml_list.set<std::string>("coarse: type", ML_CSOLVE_TYPE);
  ml_list.set<int>("coarse: max size", ML_CSOLVE_MAX_SIZE);

  return dpc_list;
}


/* ******************************************************************
* Block ILU preconditioner sublist
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateBILU_List_(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList bilu_list;

  bilu_list.set<std::string>("preconditioner type", "block ilu");

  double bilu_relax_value(ILU_RLXVAL);
  double bilu_abs_thresh(ILU_ABSTHR);
  double bilu_rel_thresh(ILU_RELTHR);
  int bilu_level_of_fill(ILU_LVLFILL);
  int bilu_overlap(ILU_OLV);
  std::string nonlinear_solver("NKA");

  if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
    Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters");
    if (ncp_list.isSublist("Unstructured Algorithm")) {
      Teuchos::ParameterList& ncpu_list = ncp_list.sublist("Unstructured Algorithm");
      if (ncpu_list.isSublist("Preconditioners")) {
        Teuchos::ParameterList& ncpup_list = ncpu_list.sublist("Preconditioners");
        if (ncpup_list.isSublist("Block ILU")) {
          Teuchos::ParameterList& ilu_list = ncpup_list.sublist("Block ILU");
          if (ilu_list.isParameter("Block ILU relax value")) {
            bilu_relax_value = ilu_list.get<double>("Block ILU relax value");
          }
          if (ilu_list.isParameter("Block ILU relative threshold")) {
            bilu_rel_thresh = ilu_list.get<double>("Block ILU relative threshold");
          }
          if (ilu_list.isParameter("Block ILU absolute threshold")) {
            bilu_abs_thresh = ilu_list.get<double>("Block ILU absolute threshold");
          }
          if (ilu_list.isParameter("Block ILU level of fill")) {
            bilu_level_of_fill = ilu_list.get<int>("Block ILU level of fill");
          }
          if (ilu_list.isParameter("Block ILU overlap")) {
            bilu_overlap = ilu_list.get<int>("Block ILU overlap");
          }
        }
      }
      Teuchos::ParameterList& solver_list = ncpu_list.sublist("Nonlinear Solver");
      if (solver_list.isParameter("Nonlinear Solver Type")) {
	nonlinear_solver = solver_list.get<std::string>("Nonlinear Solver Type");
      }
      if (nonlinear_solver == std::string("Newton")) {
	bilu_list.set<std::string>("discretization method", "fv: default");
      }
      else{
	bilu_list.set<std::string>("discretization method", "mfd: optimized for sparsity");
      }
    }
  }

  Teuchos::ParameterList& p_list = bilu_list.sublist("block ilu parameters");
  p_list.set<double>("fact: relax value",bilu_relax_value);
  p_list.set<double>("fact: absolute threshold",bilu_abs_thresh);
  p_list.set<double>("fact: relative threshold",bilu_rel_thresh);
  p_list.set<int>("fact: level-of-fill",bilu_level_of_fill);
  p_list.set<int>("overlap",bilu_overlap);
  p_list.set<std::string>("schwarz: combine mode","Add");

  return bilu_list;
}


/* ******************************************************************
* HypreBoomerAMG preconditioner sublist
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateHypreAMG_List_(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList dpc_list;
  dpc_list.set<std::string>("preconditioner type", "boomer amg");

  double tol(AMG_TOL);
  int ncycles(AMG_NCYC);
  int nsmooth(AMG_NSMOOTH);
  double strong_threshold(AMG_STR_THR);
  std::string nonlinear_solver("NKA");

  if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
    Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters");
    if (ncp_list.isSublist("Unstructured Algorithm")) {
      Teuchos::ParameterList& ncpu_list = ncp_list.sublist("Unstructured Algorithm");
      if (ncpu_list.isSublist("Preconditioners")) {
        Teuchos::ParameterList& ncpup_list = ncpu_list.sublist("Preconditioners");
        if (ncpup_list.isSublist("Hypre AMG")) {
          Teuchos::ParameterList& hypre_list = ncpup_list.sublist("Hypre AMG");
          if (hypre_list.isParameter("Hypre AMG cycle applications")) {
            ncycles = hypre_list.get<int>("Hypre AMG cycle applications");
          }
          if (hypre_list.isParameter("Hypre AMG smoother sweeps")) {
            nsmooth = hypre_list.get<int>("Hypre AMG smoother sweeps");
          }
          if (hypre_list.isParameter("Hypre AMG tolerance")) {
            tol = hypre_list.get<double>("Hypre AMG tolerance");
          }
          if (hypre_list.isParameter("Hypre AMG strong threshold")) {
            strong_threshold = hypre_list.get<double>("Hypre AMG strong threshold");
          }
        }
      }
      Teuchos::ParameterList& solver_list = ncpu_list.sublist("Nonlinear Solver");
      if (solver_list.isParameter("Nonlinear Solver Type")) {
	nonlinear_solver = solver_list.get<std::string>("Nonlinear Solver Type");
      }
      if (nonlinear_solver == std::string("Newton")) {
	dpc_list.set<std::string>("discretization method", "fv: default");
      }
      else{
	dpc_list.set<std::string>("discretization method", "mfd: optimized for sparsity");
      }
    }
  }

  Teuchos::ParameterList& amg_list = dpc_list.sublist("boomer amg parameters");
  amg_list.set<double>("tolerance", tol);
  amg_list.set<int>("smoother sweeps", nsmooth);
  amg_list.set<int>("cycle applications", ncycles);
  amg_list.set<double>("strong threshold", strong_threshold);
  amg_list.set<int>("cycle type", 1);
  amg_list.set<int>("coarsen type", 0);
  amg_list.set<int>("verbosity", 0);
  if (flow_single_phase) {
    amg_list.set<int>("relaxation type down", 3);
    amg_list.set<int>("relaxation type up", 4);
  } else {
    amg_list.set<int>("relaxation type", 3);
  }

  return dpc_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
