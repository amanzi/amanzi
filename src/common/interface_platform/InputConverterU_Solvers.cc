/*
  This is the input component of the Amanzi code. 

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

  // define defaults...
  double tol = LINEAR_SOLVER_TOL;
  int maxiter = LINEAR_SOLVER_MAXITER;
  std::string method = LINEAR_SOLVER_METHOD;
  std::string prec = LINEAR_SOLVER_PC;

  // get values from Execution control list if they exist
  DOMNodeList* node_list = doc_->getElementsByTagName(mm.transcode("unstr_linear_solver"));
  if (node_list->getLength() > 0) {
    DOMNodeList* children = node_list->item(0)->getChildNodes();

    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i);
      if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
        char* tagname = mm.transcode(inode->getNodeName());
        char* text_content = mm.transcode(inode->getTextContent());

        if (strcmp(tagname, "tolerance") == 0) {
          tol = std::strtod(text_content, NULL);
        } else if (strcmp(tagname, "max_iterations") == 0) {
          maxiter = std::strtol(text_content, NULL, 10);
        } else if (strcmp(tagname, "method") == 0) {
          method = TrimString_(text_content);
        } else if (strcmp(tagname, "preconditioner") == 0) {
          prec = TrimString_(text_content);
        }
      }
    }
  }

  // Aztec solver
  Teuchos::ParameterList& aztecoo_list = out_list.sublist("AztecOO");
  aztecoo_list.set<std::string>("preconditioner", prec);
  aztecoo_list.set<std::string>("iterative method", method);
  {
    method.append(" parameters");
    Teuchos::ParameterList& method_list = aztecoo_list.sublist(method);
    method_list.set<double>("error tolerance", tol);
    method_list.set<int>("maximum number of iterations", maxiter);
    method_list.set<int>("controller training start", 0);  // two gmres extensions
    method_list.set<int>("controller training end", 3);
    method_list.sublist("VerboseObject") = verb_list_.sublist("VerboseObject");
  }

  // add default PCG solver
  Teuchos::ParameterList& pcg_list = out_list.sublist("PCG with Hypre AMG");
  pcg_list.set<std::string>("preconditioner", prec);
  pcg_list.set<std::string>("iterative method", "pcg");
  {
    Teuchos::ParameterList& method_list = pcg_list.sublist("pcg parameters");
    method_list.set<double>("error tolerance", tol);
    method_list.set<int>("maximum number of iterations", maxiter);
    method_list.sublist("VerboseObject") = verb_list_.sublist("VerboseObject");
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
      method_list.sublist("VerboseObject").set<std::string>("Verbosity Level", "low");
    }
  }

  return out_list;
}


/* ******************************************************************
* ML preconditioner sublist
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateTrilinosML_()
{
  Teuchos::ParameterList out_list;
  out_list.set<std::string>("preconditioner type", "ml");

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
  out_list.set<std::string>("preconditioner type", "block ilu");

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
  out_list.set<std::string>("preconditioner type", "boomer amg");

  MemoryManager mm;

  // default parameters that we try to override
  double tol(HYPRE_AMG_TOL);
  int ncycles(HYPRE_AMG_NCYC);
  int nsmooth(HYPRE_AMG_NSMOOTH);
  double strong_threshold(HYPRE_AMG_STR_THR);

  bool flag;
  DOMNode* node = GetUniqueElementByTagsString_("unstr_preconditioners, hypre_amg", flag);

  if (flag) {
    DOMNodeList* children = node->getChildNodes();

    for (int i = 0; i < children->getLength(); i++) {
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
  amg_list.set<int>("verbosity", 0);
  if (flow_single_phase_) {
    amg_list.set<int>("relaxation type down", 3);
    amg_list.set<int>("relaxation type up", 4);
  } else {
    amg_list.set<int>("relaxation type", 3);
  }

  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
