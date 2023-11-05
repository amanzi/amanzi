/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Input Converter

*/

#include <algorithm>
#include <sstream>
#include <string>

//TPLs
#include <xercesc/dom/DOM.hpp>

// Amanzi's
#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"
#include "StringExt.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Create operators sublist
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateTimeIntegrator_(const std::string& err_options,
                                          const std::string& nonlinear_solver,
                                          bool modify_correction,
                                          const std::string& controls,
                                          const std::string& linsolver,
                                          double dt_cut_default,
                                          double dt_inc_default)
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating time integrator" << std::endl;

  MemoryManager mm;
  DOMNode* node;

  // error control options
  std::vector<std::string> tmp = CharToStrings_(err_options.c_str());
  out_list.set<Teuchos::Array<std::string>>("error control options", tmp);

  // linear solver
  bool flag;
  std::string prec(TI_PRECONDITIONER);
  node = GetUniqueElementByTagsString_(controls + ", preconditioner", flag);
  if (flag) prec = mm.transcode(node->getTextContent());

  out_list.set<std::string>("linear solver", linsolver);
  out_list.set<std::string>("preconditioner", prec);
  // out_list.set<std::string>("preconditioner enhancement", "none");

  // pressure-lambda constraints
  Teuchos::ParameterList& plamb = out_list.sublist("pressure-lambda constraints");
  plamb.set<std::string>("method", "projection");
  plamb.set<bool>("inflow krel correction", true);
  plamb.set<std::string>("linear solver", TI_PLAMBDA_SOLVER);

  // time stepping method for high-level PK only
  out_list.set<std::string>("time integration method", "BDF1");
  Teuchos::ParameterList& bdf1 = out_list.sublist("BDF1");

  // use standard timestep controller type
  std::string name(TI_TIMESTEP_CONTROLLER);
  node = GetUniqueElementByTagsString_(controls + ", timestep_controller", flag);
  if (flag) name = mm.transcode(node->getTextContent());

  bdf1.set<std::string>("timestep controller type", name);
  Teuchos::ParameterList& controller = bdf1.sublist("timestep controller " + name + " parameters");
  controller.set<int>("max iterations", TI_MAX_ITERATIONS)
    .set<int>("min iterations", TI_MIN_ITERATIONS)
    .set<double>("time step increase factor", dt_inc_default)
    .set<double>("time step reduction factor", dt_cut_default)
    .set<double>("max time step", MAXIMUM_TIMESTEP)
    .set<double>("min time step", MINIMUM_TIMESTEP);
  if (name == "adaptive")
    controller.set<double>("relative tolerance", 1e-4).set<double>("absolute tolerance", 10.0);

  // nonlinear solver
  Teuchos::ParameterList* solver;

  if (nonlinear_solver == std::string("newton") ||
      nonlinear_solver == std::string("newton-picard")) {
    bdf1.set<std::string>("solver type", "Newton");
    Teuchos::ParameterList& test = bdf1.sublist("Newton parameters");
    solver = &test;
    solver->set<double>("nonlinear tolerance", NONLINEAR_TOLERANCE);
    solver->set<double>("diverged tolerance", NKA_DIVERG_TOL);
    solver->set<double>("max du growth factor", INC_DIVERG_FACTOR);
    solver->set<int>("max divergent iterations", MAX_DIVERG_ITERATIONS);
    solver->set<int>("limit iterations", NKA_LIMIT_ITERATIONS);
    solver->set<bool>("modify correction", modify_correction);
    solver->set<std::string>("monitor", "monitor update");
  } else if (nonlinear_solver == "nka") {
    bdf1.set<std::string>("solver type", "nka");
    Teuchos::ParameterList& test = bdf1.sublist("nka parameters");
    solver = &test;
    solver->set<double>("nonlinear tolerance", NONLINEAR_TOLERANCE);
    solver->set<double>("diverged tolerance", NKA_DIVERG_TOL);
    solver->set<double>("diverged l2 tolerance", NKA_DIVERG_TOL);
    solver->set<double>("max du growth factor", INC_DIVERG_FACTOR);
    solver->set<int>("max divergent iterations", MAX_DIVERG_ITERATIONS);
    solver->set<int>("max nka vectors", NKA_NUM_VECTORS);
    solver->set<int>("limit iterations", NKA_LIMIT_ITERATIONS);
    solver->set<bool>("modify correction", modify_correction);
  } else if (nonlinear_solver == std::string("jfnk")) {
    bdf1.set<std::string>("solver type", "JFNK");
    Teuchos::ParameterList& test = bdf1.sublist("JFNK parameters");
    solver = &test;
    solver->set<double>("typical solution value", 1.0);

    Teuchos::ParameterList& list_tmp = solver->sublist("nonlinear solver");
    list_tmp.set<std::string>("solver type", "Newton");
    Teuchos::ParameterList& newton = list_tmp.sublist("Newton parameters");
    newton.set<double>("diverged tolerance", NKA_DIVERG_TOL);
    newton.set<double>("max du growth factor", INC_DIVERG_FACTOR);
    newton.set<int>("max divergent iterations", MAX_DIVERG_ITERATIONS);
    newton.set<int>("max nka vectors", NKA_NUM_VECTORS);
    newton.set<int>("limit iterations", NKA_LIMIT_ITERATIONS);

    Teuchos::ParameterList& jfmat = solver->sublist("JF matrix parameters");
    jfmat.set<double>("finite difference epsilon", 1.0e-8);
    jfmat.set<std::string>("method for epsilon", "Knoll-Keyes");

    Teuchos::ParameterList& linop = solver->sublist("linear operator");
    linop.set<std::string>("iterative method", "gmres");
    Teuchos::ParameterList& gmres = linop.sublist("gmres parameters");
    gmres.set<double>("error tolerance", 1e-7);
    gmres.set<int>("maximum number of iterations", 100);
    std::vector<std::string> criteria;
    criteria.push_back("relative rhs");
    criteria.push_back("relative residual");
    gmres.set<Teuchos::Array<std::string>>("convergence criteria", criteria);
  } else {
    Errors::Message msg;
    msg << "In the definition of \"unstr_nonlinear_solver\" you must specify either"
        << " 'nka', 'newton', 'jfnk', or 'newton-picard'.\n";
    Exceptions::amanzi_throw(msg);
  }

  // remaining BDF1 parameters
  bdf1.set<int>("max preconditioner lag iterations", TI_MAX_PC_LAG);
  bdf1.set<bool>("extrapolate initial guess", true);
  bdf1.set<double>("restart tolerance relaxation factor", TI_TOL_RELAX_FACTOR);
  bdf1.set<double>("restart tolerance relaxation factor damping", TI_TOL_RELAX_FACTOR_DAMPING);
  bdf1.set<int>("nonlinear iteration initial guess extrapolation order", 1);

  // controller options
  {
    std::vector<std::string> options({ "max_iterations", "min_iterations" });
    for (auto opt : options) {
      node = GetUniqueElementByTagsString_(controls + ", " + opt, flag);
      if (flag) {
        auto tmp(opt);
        std::replace(tmp.begin(), tmp.end(), '_', ' ');
        controller.set<int>(tmp, strtol(mm.transcode(node->getTextContent()), NULL, 10));
      }
    }
  }

  {
    std::vector<std::string> options({ "time_step_reduction_factor", "time_step_increase_factor" });
    for (auto opt : options) {
      node = GetUniqueElementByTagsString_(controls + ", " + opt, flag);
      if (flag) {
        auto tmp(opt);
        std::replace(tmp.begin(), tmp.end(), '_', ' ');
        controller.set<double>(tmp, strtod(mm.transcode(node->getTextContent()), NULL));
      }
    }
  }

  // solver options
  {
    std::vector<std::string> options({ "limit_iterations", "max_divergent_iterations" });
    for (auto opt : options) {
      node = GetUniqueElementByTagsString_(controls + ", " + opt, flag);
      if (flag) {
        auto tmp(opt);
        std::replace(tmp.begin(), tmp.end(), '_', ' ');
        solver->set<int>(tmp, strtol(mm.transcode(node->getTextContent()), NULL, 10));
      }
    }
  }

  {
    std::vector<std::string> options(
      { "nonlinear_tolerance", "diverged_tolerance", "max_error_growth_factor" });
    for (auto opt : options) {
      node = GetUniqueElementByTagsString_(controls + ", " + opt, flag);
      if (flag) {
        auto tmp(opt);
        std::replace(tmp.begin(), tmp.end(), '_', ' ');
        solver->set<double>(tmp, strtod(mm.transcode(node->getTextContent()), NULL));
      }
    }
  }

  node = GetUniqueElementByTagsString_(controls + ", nonlinear_iteration_divergence_factor", flag);
  if (flag)
    solver->set<double>("max du growth factor", strtod(mm.transcode(node->getTextContent()), NULL));

  node = GetUniqueElementByTagsString_(controls + ", monitor", flag);
  if (flag) solver->set<std::string>("monitor", mm.transcode(node->getTextContent()));


  // bdf1 options
  {
    std::vector<std::string> options({ "max_preconditioner_lag_iterations",
                                       "nonlinear_iteration_initial_guess_extrapolation_order" });
    for (auto opt : options) {
      node = GetUniqueElementByTagsString_(controls + ", " + opt, flag);
      if (flag) {
        auto tmp(opt);
        std::replace(tmp.begin(), tmp.end(), '_', ' ');
        bdf1.set<int>(tmp, strtol(mm.transcode(node->getTextContent()), NULL, 10));
      }
    }
  }

  {
    std::vector<std::string> options(
      { "restart_tolerance_relaxation_factor", "restart_tolerance_relaxation_factor_damping" });
    for (auto opt : options) {
      node = GetUniqueElementByTagsString_(controls + ", " + opt, flag);
      if (flag) {
        auto tmp(opt);
        std::replace(tmp.begin(), tmp.end(), '_', ' ');
        bdf1.set<double>(tmp, strtod(mm.transcode(node->getTextContent()), NULL));
      }
    }
  }

  // other options
  node = GetUniqueElementByTagsString_(controls + ", error_control_options", flag);
  if (flag)
    out_list.set<Teuchos::Array<std::string>>("error control options",
                                              CharToStrings_(mm.transcode(node->getTextContent())));

  node = GetUniqueElementByTagsString_(controls + ", preconditioner", flag);
  if (flag) {
    std::string text = GetTextContentS_(node, "hypre_amg, trilinos_ml, block_ilu");
    if (text == "hypre_amg") text = "Hypre AMG";
    if (text == "trilinos_ml") text = "Trilinos ML";
    if (text == "block_ilu") text = "Block ILU";
    out_list.set<std::string>("preconditioner", text);
  }

  // initialization
  node = GetUniqueElementByTagsString_(controls + ", unstr_initialization", flag);
  if (flag) {
    Teuchos::ParameterList& init = out_list.sublist("initialization");
    init = TranslateInitialization_(controls);
  }

  // overwrite parameters for special solvers
  if (nonlinear_solver == "newton" || nonlinear_solver == "newton-picard") {
    std::stringstream ss;
    ss << "GMRES for Newton-" << gmres_solvers_.size();

    bdf1.set<int>("max preconditioner lag iterations", 0);
    bdf1.set<bool>("extrapolate initial guess", false);
    out_list.set<std::string>("linear solver", ss.str());
    out_list.set<std::string>("preconditioner enhancement", ss.str());

    double nonlinear_tol = solver->get<double>("nonlinear tolerance");
    gmres_solvers_.push_back(std::make_pair(ss.str(), nonlinear_tol));
  }

  bdf1.sublist("verbose object") = verb_list_.sublist("verbose object");
  return out_list;
}


/* ******************************************************************
* Translate initialization sublist for time integrator
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateInitialization_(const std::string& unstr_controls)
{
  Teuchos::ParameterList out_list;

  MemoryManager mm;
  DOMNode* node;

  // set defaults
  out_list.set<std::string>("method", "saturated solver");
  out_list.set<std::string>("linear solver", TI_SOLVER);

  // overwite defaults using numerical controls
  bool flag;
  std::string method;
  std::string controls(unstr_controls + ", unstr_initialization");

  node = GetUniqueElementByTagsString_(controls + ", method", flag);
  if (flag) {
    method = GetTextContentS_(node, "picard, darcy_solver");
    if (method == "darcy_solver") method = "saturated solver";
    out_list.set<std::string>("method", method);
  }

  node = GetUniqueElementByTagsString_(controls + ", wells_status", flag);
  if (flag) {
    std::string text = mm.transcode(node->getTextContent());
    flag = (text == "on");
  }
  out_list.set<bool>("active wells", flag);

  node = GetUniqueElementByTagsString_(controls + ", clipping_saturation", flag);
  if (flag)
    out_list.set<double>("clipping saturation value",
                         strtod(mm.transcode(node->getTextContent()), NULL));

  node = GetUniqueElementByTagsString_(controls + ", clipping_pressure", flag);
  if (flag)
    out_list.set<double>("clipping pressure value",
                         strtod(mm.transcode(node->getTextContent()), NULL));

  node = GetUniqueElementByTagsString_(controls + ", linear_solver", flag);
  if (flag) {
    std::string text = mm.transcode(node->getTextContent());
    if (text == "aztec00" || text == "aztecoo") text = "AztecOO";
    out_list.set<std::string>("linear solver", text);
  }

  if (method == "picard") {
    Teuchos::ParameterList& pic_list = out_list.sublist("picard parameters");
    out_list.set<std::string>("linear solver", "GMRES with Hypre AMG");
    pic_list.set<double>("convergence tolerance", PICARD_TOLERANCE);
    pic_list.set<int>("maximum number of iterations", PICARD_MAX_ITERATIONS);

    node = GetUniqueElementByTagsString_(controls + ", convergence_tolerance", flag);
    if (flag)
      pic_list.set<double>("convergence tolerance",
                           strtod(mm.transcode(node->getTextContent()), NULL));

    node = GetUniqueElementByTagsString_(controls + ", max_iterations", flag);
    if (flag)
      pic_list.set<int>("maximum number of iterations",
                        strtol(mm.transcode(node->getTextContent()), NULL, 10));
  }

  return out_list;
}


/* ******************************************************************
* Create operators sublist:
*   disc_methods     = {primary, secondary}
*   pc_method        = "linerized_operator" | "diffusion_operator"
*   nonlinear_solver = "" | "Newton"
*   nonlinear_coef   = "upwind-darcy_velocity" | etc
*   extentions       = "" | "vapor matrix"
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateDiffusionOperator_(const std::string& disc_methods,
                                             const std::string& pc_method,
                                             const std::string& nonlinear_solver,
                                             const std::string& nonlinear_coef,
                                             const std::string& extensions,
                                             const std::string& domain,
                                             bool gravity,
                                             const std::string& pk)
{
  Teuchos::ParameterList out_list;
  Teuchos::ParameterList tmp_list;

  // process primary and secondary discretization methods
  std::vector<std::string> methods = CharToStrings_(disc_methods.c_str());
  for (int i = 0; i < methods.size(); ++i) {
    std::string tmp(methods[i]);
    Amanzi::replace_all(tmp, "-", ": ");
    replace(tmp.begin(), tmp.end(), '_', ' ');
    if (tmp == "mfd: two point flux approximation") tmp = "mfd: two-point flux approximation";
    methods[i] = tmp;
  }
  if (methods.size() == 1) methods.push_back("mfd: optimized for sparsity");

  tmp_list.set<std::string>("discretization primary", methods[0]);
  tmp_list.set<std::string>("discretization secondary", methods[1]);
  if (gravity_on_) tmp_list.set<bool>("gravity", gravity);

  // process nonlinear coefficient for PDE operators
  std::string nonlinear_coef_out(nonlinear_coef);

  Amanzi::replace_all(nonlinear_coef_out, "-", ": ");
  std::replace(nonlinear_coef_out.begin(), nonlinear_coef_out.end(), '_', ' ');

  if (nonlinear_coef == "upwind-darcy_velocity")
    nonlinear_coef_out = "upwind: face";
  else if (nonlinear_coef == "upwind-amanzi" || nonlinear_coef == "upwind-amanzi_new")
    nonlinear_coef_out = "divk: cell-face";

  // -- limitations
  if (fracture_network_ && domain == "fracture" && pk != "multiphase")
    nonlinear_coef_out = "standard: cell";

  if (nonlinear_coef_out != "") tmp_list.set("nonlinear coefficient", nonlinear_coef_out);

  // process schema
  if (methods[0] != "fv: default" && methods[0] != "nlfv: default") {
    Teuchos::Array<std::string> stensil(2);
    stensil[0] = "face";
    stensil[1] = "cell";
    tmp_list.set<Teuchos::Array<std::string>>("schema", stensil);

    if (pc_method != "linearized_operator" && fracture_regions_.size() == 0) stensil.remove(1);
    tmp_list.set<Teuchos::Array<std::string>>("preconditioner schema", stensil);
    if (gravity && nonlinear_coef == "upwind-amanzi_new")
      tmp_list.set<std::string>("gravity term discretization", "finite volume");
  } else {
    Teuchos::Array<std::string> stensil(1);
    stensil[0] = "cell";
    tmp_list.set<Teuchos::Array<std::string>>("schema", stensil);
    tmp_list.set<Teuchos::Array<std::string>>("preconditioner schema", stensil);
  }

  // fractured matrix
  if (fracture_network_ && domain != "fracture")
    tmp_list.set<Teuchos::Array<std::string>>("fracture", fracture_regions_);

  // create two operators for matrix and preconditioner
  // Note that PK may use only one of them.
  out_list.sublist("diffusion operator").sublist("matrix") = tmp_list;
  out_list.sublist("diffusion operator").sublist("preconditioner") = tmp_list;

  // extensions
  if (extensions == "vapor matrix") {
    Teuchos::ParameterList& vapor = out_list.sublist("diffusion operator").sublist("vapor matrix");
    vapor = tmp_list;
    vapor.set<std::string>("nonlinear coefficient", "standard: cell");
    vapor.set<bool>("exclude primary terms", false);
    vapor.set<bool>("scaled constraint equation", false);
    if (gravity_on_) vapor.set<bool>("gravity", "false");
    vapor.set<std::string>("Newton correction", "none");
    if (fracture_network_ && domain != "fracture")
      vapor.set<Teuchos::Array<std::string>>("fracture", fracture_regions_);
  }

  // fixing miscalleneous scenarious
  if (pc_method == "linearized_operator" && pk == "flow") {
    out_list.sublist("diffusion operator")
      .sublist("preconditioner")
      .set<std::string>("Newton correction", "approximate Jacobian");
  }

  if (nonlinear_solver == "newton") {
    Teuchos::ParameterList& pc_list =
      out_list.sublist("diffusion operator").sublist("preconditioner");
    pc_list.set<std::string>("Newton correction", "true Jacobian");
  } else if (nonlinear_solver == "newton-picard") {
    Teuchos::ParameterList& pc_list =
      out_list.sublist("diffusion operator").sublist("preconditioner");
    pc_list.set<std::string>("Newton correction", "approximate Jacobian");
  }

  return out_list;
}

} // namespace AmanziInput
} // namespace Amanzi
