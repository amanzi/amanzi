/*
  Input Converter

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <sstream>
#include <string>

//TPLs
#include <xercesc/dom/DOM.hpp>

// Amanzi's
#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Collects default preconditioners
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslatePreconditioners_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Translating preconditioners" << std::endl;
  }

  out_list.sublist("Trilinos ML") = TranslateTrilinosML_();
  out_list.sublist("Hypre AMG") = TranslateHypreAMG_();
  out_list.sublist("Block ILU") = TranslateBILU_();
  if (multiphase_) {
    out_list.sublist("Euclid") = TranslateEuclid_();
  }

  return out_list;
}


/* ******************************************************************
* Collects linear solvers
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateSolvers_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Translating solvers" << std::endl;
  }

  MemoryManager mm;

  // User's solver
  out_list.sublist("AztecOO") = TranslateLinearSolvers_("", LINEAR_SOLVER_METHOD, "");

  // add PCG and GMRES solvers (generic or specialized)
  out_list.sublist("Dispersion Solver") = TranslateLinearSolvers_(
      "unstr_transport_controls, dispersion_linear_solver", "pcg", "");
  out_list.sublist("Dispersion Solver").sublist("pcg parameters")
          .sublist("verbose object").set<std::string>("verbosity level", "low");

  out_list.sublist("PCG with Hypre AMG") = TranslateLinearSolvers_(
      "unstr_flow_controls, saturated_linear_solver", "pcg", "");

  if (pk_model_.find("flow") != pk_model_.end()) {
    std::string enforce = (pk_model_["flow"] == "richards") ? "gmres" : "";
    out_list.sublist("GMRES with Hypre AMG") = TranslateLinearSolvers_(
        "unstr_flow_controls, constraints_linear_solver", LINEAR_SOLVER_METHOD, enforce);
  }

  // add default "GMRES for Newton" solver
  for (int i = 0; i < gmres_solvers_.size(); ++i) {
    Teuchos::ParameterList& gmres_list = out_list.sublist(gmres_solvers_[i].first);
    gmres_list.set<std::string>("iterative method", "gmres");
    {
      Teuchos::ParameterList& method_list = gmres_list.sublist("gmres parameters");
      method_list.set<double>("error tolerance", gmres_solvers_[i].second * 1e-2);
      method_list.set<int>("maximum number of iterations", 50);
      std::vector<std::string> criteria;
      criteria.push_back("relative rhs");
      criteria.push_back("relative residual");
      method_list.set<Teuchos::Array<std::string> >("convergence criteria", criteria);
      method_list.set<int>("controller training start", 0);
      method_list.set<int>("controller training end", 3);
      method_list.sublist("verbose object").set<std::string>("verbosity level", "low");
      method_list.set<bool>("release Krylov vectors", true);
    }
  }

  // direct solver
  Teuchos::ParameterList& amesos_list = out_list.sublist("AMESOS");
  amesos_list.set<std::string>("direct method", "amesos");
  amesos_list.sublist("amesos parameters")
      .template set<std::string>("solver name", "basker")
      .template set<int>("amesos version", 2);

  return out_list;
}


/* ******************************************************************
* Creates linear solver list from input paramaters.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateLinearSolvers_(
    std::string tags, std::string method_default, std::string method_enforce) 
{
  Errors::Message msg;
  Teuchos::ParameterList plist;

  DOMNode* node;
  MemoryManager mm;

  int maxiter(LINEAR_SOLVER_MAXITER);
  double tol(LINEAR_SOLVER_TOL);
  std::string tags_default("unstructured_controls, unstr_linear_solver");
  std::string method(method_default), prec(LINEAR_SOLVER_PC);

  // verify that method is admissible
  bool flag;
  node = GetUniqueElementByTagsString_(tags_default + ", method", flag);
  if (flag) method = GetTextContentS_(node, "pcg, gmres"); 
  node = GetUniqueElementByTagsString_(tags + ", method", flag);
  if (flag) method = GetTextContentS_(node, "pcg, gmres"); 

  if (method_enforce != "" && method != method_enforce) {
    msg << "Expect method=\"" << method_enforce 
        << "\" under default node=\"" << tags_default 
        << "\"\n or under PK node \"" << tags << "\".\n";
    Exceptions::amanzi_throw(msg);
  }

  // collect other parameters
  node = GetUniqueElementByTagsString_(tags_default + ", tolerance", flag);
  if (flag) tol = GetTextContentD_(node, "-");
  node = GetUniqueElementByTagsString_(tags + ", tolerance", flag);
  if (flag) tol = GetTextContentD_(node, "-");

  node = GetUniqueElementByTagsString_(tags_default + ", max_iterations", flag);
  if (flag) maxiter = std::strtol(mm.transcode(node->getTextContent()), NULL, 10);
  node = GetUniqueElementByTagsString_(tags + ", max_iterations", flag);
  if (flag) maxiter = std::strtol(mm.transcode(node->getTextContent()), NULL, 10);

  // populate parameter list
  plist.set<std::string>("iterative method", method);
 
  Teuchos::ParameterList& slist = plist.sublist(method + " parameters");
  slist.set<double>("error tolerance", tol);
  slist.set<int>("maximum number of iterations", maxiter);

  // parameters without 2.x support
  if (method == "gmres") {
    slist.set<int>("controller training start", 0);
    slist.set<int>("controller training end", 3);
  }

  slist.sublist("verbose object") = verb_list_.sublist("verbose object");

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << method << ": itrs=" << maxiter 
               << ", tol=" << tol << ", pc=" << prec << std::endl;

  return plist;
}


/* ******************************************************************
* ML preconditioner sublist
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateTrilinosML_()
{
  Teuchos::ParameterList out_list;
  out_list.set<std::string>("preconditioning method", "ml");

  MemoryManager mm;

  // default parameters that we try to override
  double aggthr(TRILINOS_ML_AGG_THR);
  std::string smthtyp(TRILINOS_ML_SMOOTHER);
  int ncycles(TRILINOS_ML_NCYC);
  int nsmooth(TRILINOS_ML_NSMOOTH);

  bool flag;
  DOMNode* node = GetUniqueElementByTagsString_("unstr_preconditioners, trilinos_ml", flag);
  if (flag) {
    DOMNodeList* children = node->getChildNodes();

    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i);
      char* tagname = mm.transcode(inode->getNodeName());
      char* text_content = mm.transcode(inode->getTextContent());

      if (strcmp(tagname, "trilinos_smoother_type") == 0) {
        smthtyp = TrimString_(text_content);
      } else if (strcmp(tagname, "trilinos_threshold") == 0) {
        aggthr = std::strtod(text_content, NULL);
      } else if (strcmp(tagname, "trilinos_smoother_sweeps") == 0) {
        nsmooth = std::strtol(text_content, NULL, 10);
      } else if (strcmp(tagname, "trilinos_cycle_applications") == 0) {
        ncycles = std::strtol(text_content, NULL, 10);
      }
    }
  }

  Teuchos::ParameterList& ml_list = out_list.sublist("ml parameters");
  ml_list.set<int>("ML output", TRILINOS_ML_OUTPUT);
  ml_list.set<int>("max levels", TRILINOS_ML_MAXLVLS);
  ml_list.set<std::string>("prec type",TRILINOS_ML_PRECTYPE);
  ml_list.set<int>("cycle applications", ncycles);
  ml_list.set<std::string>("aggregation: type", TRILINOS_ML_AGGTYPE);
  ml_list.set<double>("aggregation: damping factor", TRILINOS_ML_AGGDAMP);
  ml_list.set<double>("aggregation: threshold", aggthr);
  ml_list.set<std::string>("eigen-analysis: type",TRILINOS_ML_EIGENANAL_TYPE);
  ml_list.set<int>("eigen-analysis: iterations", TRILINOS_ML_EIGENANAL_ITERS);
  ml_list.set<int>("smoother: sweeps", nsmooth);
  ml_list.set<double>("smoother: damping factor", TRILINOS_ML_SMOOTH_DAMP);
  ml_list.set<std::string>("smoother: pre or post", TRILINOS_ML_SMOOTH_PRE_POST);
  ml_list.set<std::string>("smoother: type", smthtyp);
  ml_list.set<double>("smoother: damping factor", TRILINOS_ML_SMOOTH_DAMP);
  ml_list.set<std::string>("coarse: type", TRILINOS_ML_CSOLVE_TYPE);
  ml_list.set<int>("coarse: max size", TRILINOS_ML_CSOLVE_MAX_SIZE);

  return out_list;
}


/* ******************************************************************
* Block ILU preconditioner sublist
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateBILU_()
{
  Teuchos::ParameterList out_list;
  out_list.set<std::string>("preconditioning method", "block ilu");

  MemoryManager mm;

  // default parameters that we try to override
  double bilu_relax_value(TRILINOS_ILU_RLXVAL);
  double bilu_abs_thresh(TRILINOS_ILU_ABSTHR);
  double bilu_rel_thresh(TRILINOS_ILU_RELTHR);
  int bilu_level_of_fill(TRILINOS_ILU_LVLFILL);
  int bilu_overlap(TRILINOS_ILU_OLV);

  bool flag;
  DOMNode* node = GetUniqueElementByTagsString_("unstr_preconditioners, block_ilu", flag);

  if (flag) {
    DOMNodeList* children = node->getChildNodes();

    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i);
      char* tagname = mm.transcode(inode->getNodeName());
      char* text_content = mm.transcode(inode->getTextContent());

      if (strcmp(tagname, "ilu_overlap") == 0) {
        bilu_overlap = std::strtol(text_content, NULL, 10);
      } else if (strcmp(tagname, "ilu_relax") == 0) {
        bilu_relax_value = std::strtod(text_content, NULL);
      } else if (strcmp(tagname, "ilu_rel_threshold") == 0) {
        bilu_rel_thresh = std::strtod(text_content, NULL);
      } else if (strcmp(tagname, "ilu_abs_threshold") == 0) {
        bilu_abs_thresh = std::strtod(text_content, NULL);
      } else if (strcmp(tagname, "ilu_level_of_fill") == 0) {
        bilu_level_of_fill = std::strtol(text_content, NULL, 10);
      } 
    }
  }

  Teuchos::ParameterList& p_list = out_list.sublist("block ilu parameters");
  p_list.set<double>("fact: relax value", bilu_relax_value);
  p_list.set<double>("fact: absolute threshold", bilu_abs_thresh);
  p_list.set<double>("fact: relative threshold", bilu_rel_thresh);
  p_list.set<int>("fact: level-of-fill", bilu_level_of_fill);
  p_list.set<int>("overlap", bilu_overlap);
  p_list.set<std::string>("schwarz: combine mode", "Add");

  return out_list;
}


/* ******************************************************************
* HypreBoomerAMG preconditioner sublist
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateHypreAMG_()
{
  Teuchos::ParameterList out_list;
  out_list.set<std::string>("preconditioning method", "boomer amg");

  MemoryManager mm;

  // default parameters that we try to override
  double tol(HYPRE_AMG_TOL);
  int ncycles(HYPRE_AMG_NCYC);
  int nsmooth(HYPRE_AMG_NSMOOTH);
  double strong_threshold(HYPRE_AMG_STR_THR);

  bool flag, block_indices(false);
  DOMNode* node = GetUniqueElementByTagsString_("unstr_preconditioners, hypre_amg", flag);

  if (flag) {
    DOMNodeList* children = node->getChildNodes();
    int nchildren = children->getLength();

    for (int i = 0; i < nchildren; i++) {
      DOMNode* inode = children->item(i);
      char* tagname = mm.transcode(inode->getNodeName());
      char* text_content = mm.transcode(inode->getTextContent());

      if (strcmp(tagname, "hypre_cycle_applications") == 0) {
        ncycles = std::strtol(text_content, NULL, 10);
      } else if (strcmp(tagname, "hypre_smoother_sweeps") == 0) {
        nsmooth = std::strtol(text_content, NULL, 10);
      } else if (strcmp(tagname, "hypre_tolerance") == 0) {
        tol = std::strtod(text_content, NULL);
      } else if (strcmp(tagname, "hypre_strong_threshold") == 0) {
        strong_threshold = std::strtod(text_content, NULL);
      } else if (strcmp(tagname, "use_block_indices") == 0) {
        block_indices = (strcmp(text_content, "true") == 0);
      }
    }
  }

  Teuchos::ParameterList& amg_list = out_list.sublist("boomer amg parameters");
  amg_list.set<double>("tolerance", tol);
  amg_list.set<int>("smoother sweeps", nsmooth);
  amg_list.set<int>("cycle applications", ncycles);
  amg_list.set<double>("strong threshold", strong_threshold);
  amg_list.set<int>("cycle type", 1);
  amg_list.set<int>("coarsen type", 0);
  if (block_indices)
    amg_list.set<bool>("use block indices", block_indices);
  amg_list.set<int>("verbosity", 0);
  if (flow_single_phase_) {
    amg_list.set<int>("relaxation type down", 13);
    amg_list.set<int>("relaxation type up", 14);
  } else {
    amg_list.set<int>("relaxation type", 3);
  }

  return out_list;
}


/* ******************************************************************
* Euclid sublist
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateEuclid_()
{
  Teuchos::ParameterList out_list;
  out_list.set<std::string>("preconditioning method", "euclid");

  out_list.sublist("euclid parameters")
      .set<int>("ilu(k) fill level", 5)
      // .set<double>("ILUT drop tolerance", 0.000001)
      .set<bool>("rescale rows", true)
      .set<int>("verbosity", 0);

  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
