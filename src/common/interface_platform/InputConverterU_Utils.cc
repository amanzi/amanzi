/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (original version)
           Konstantin Lipnikov (lipnikov@lanl.gov)
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
#include "Units.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Translate list of boundary conditions.
****************************************************************** */
BCs
InputConverterU::ParseCondList_(DOMNode* node,
                                double vmin,
                                double vmax,
                                const std::string& unit,
                                bool is_bc)
{
  bool flag;
  std::string bctype;

  std::vector<DOMNode*> same_list = GetSameChildNodes_(node, bctype, flag, true);
  return ParseCondList_(same_list, bctype, vmin, vmax, unit, is_bc);
}


/* ******************************************************************
* Translate list of boundary conditions.
****************************************************************** */
BCs
InputConverterU::ParseCondList_(std::vector<DOMNode*> same_list,
                                const std::string& bctype,
                                double vmin,
                                double vmax,
                                const std::string& unit,
                                bool is_bc)
{
  int nlist = same_list.size();
  std::map<double, double> tp_values, tp_fluxes;
  std::map<double, std::string> tp_forms, tp_formulas;
  DOMElement* element;
  BCs bcs;

  for (int j = 0; j < nlist; ++j) {
    element = static_cast<DOMElement*>(same_list[j]);

    // allowed attributes
    if (HasAttribute_(element, "filename") && nlist == 1) {
      bcs.filename = GetAttributeValueS_(element, "filename", TYPE_NONE, false);
      bcs.xheader = GetAttributeValueS_(element, "times", TYPE_NONE, false, "");
      bcs.yheader = GetAttributeValueS_(element, "values", TYPE_NONE, false, "");
      bcs.variable = GetAttributeValueS_(element, "variable", TYPE_NONE, false, "");
    } else {
      double t0 = GetAttributeValueD_(element, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s", is_bc);
      tp_values[t0] = 0.0;
      tp_fluxes[t0] = 0.0;
      tp_forms[t0] = "";

      if (HasAttribute_(element, "formula")) {
        tp_formulas[t0] = GetAttributeValueS_(element, "formula", TYPE_NONE, false, "");
        tp_forms[t0] = "general";
      } else if (HasAttribute_(element, "value")) {
        tp_values[t0] =
          GetAttributeValueD_(element, "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, unit, false);
        if (is_bc) 
          tp_forms[t0] = GetAttributeValueS_(same_list[j], "function", "constant, linear, uniform");
      } else if (HasAttribute_(element, "inward_mass_flux")) {
        tp_fluxes[t0] = GetAttributeValueD_(
          element, "inward_mass_flux", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, unit, false);
      } else {
        Errors::Message msg;
        msg << "Unknown or ill-formed boundary conditions.\"\n";
        Exceptions::amanzi_throw(msg);
      }
    }
  }

  // create vectors of values and forms
  for (auto it = tp_values.begin(); it != tp_values.end(); ++it) {
    bcs.times.push_back(it->first);
    bcs.values.push_back(it->second);
    bcs.fluxes.push_back(tp_fluxes[it->first]);
    bcs.forms.push_back(tp_forms[it->first]);
    bcs.formulas.push_back(tp_formulas[it->first]);
  }

  bcs.type = bctype;
  return bcs;
}


/* ******************************************************************
* Translate function Gaussian.
* Data format: d0 exp(-|x-d1|^2 / (2 d4^2)) where d1 is space vector.
****************************************************************** */
void
InputConverterU::TranslateFunctionGaussian_(const std::vector<double>& data,
                                            Teuchos::ParameterList& bcfn)
{
  if (data.size() != dim_ + 2) {
    Errors::Message msg;
    msg << "Gaussian function requires " << dim_ + 2 << " parameters.\n";
    Exceptions::amanzi_throw(msg);
  }

  std::vector<double> data_tmp(data);

  bcfn.sublist("function-composition")
    .sublist("function1")
    .sublist("function-standard-math")
    .set<std::string>("operator", "exp")
    .set<double>("amplitude", data_tmp[0]);

  double sigma = data_tmp[dim_ + 1];
  double factor = -0.5 / sigma / sigma;
  data_tmp.pop_back();

  Teuchos::ParameterList& bc_tmp =
    bcfn.sublist("function-composition").sublist("function2").sublist("function-multiplicative");

  std::vector<double> metric(data_tmp.size(), 1.0);
  metric[0] = 0.0; // ignore time distance
  data_tmp[0] = 0.0;
  bc_tmp.sublist("function1")
    .sublist("function-squaredistance")
    .set<Teuchos::Array<double>>("x0", data_tmp)
    .set<Teuchos::Array<double>>("metric", metric);
  bc_tmp.sublist("function2").sublist("function-constant").set<double>("value", factor);
}


/* ******************************************************************
* Translate linear function.
****************************************************************** */
void InputConverterU::TranslateFunctionGradient_(double refv, 
                                                 std::vector<double>& grad,
                                                 std::vector<double>& refc,
                                                 Teuchos::ParameterList& bcfn)
{
  grad.insert(grad.begin(), 0.0);
  refc.insert(refc.begin(), 0.0);

  bcfn.sublist("function-linear")
    .set<double>("y0", refv)
    .set<Teuchos::Array<double>>("x0", refc)
    .set<Teuchos::Array<double>>("gradient", grad);
}


/* ******************************************************************
* Generic output of math data.
****************************************************************** */
bool
InputConverterU::TranslateGenericMath_(const BCs& bcs,
                                       Teuchos::ParameterList& bcfn)
{
  bool flag(false);

  if (bcs.times.size() == 1) {
    if (bcs.forms[0] == "general") {
      bcfn.sublist("function-exprtk")
        .set<int>("number of arguments", dim_ + 1)
        .set<std::string>("formula", bcs.formulas[0]);
      flag = true;
    } else {
      bcfn.sublist("function-constant").set<double>("value", bcs.values[0]);
      flag = true;
    }
  } else {
    auto forms_copy = bcs.forms;
    forms_copy.pop_back();

    for (int i = 0; i < forms_copy.size(); ++i) {
      if (forms_copy[i] == "general") {
        std::string user_fnc = "USER_FNC" + std::to_string(i);
        bcfn.sublist("function-tabular")
          .sublist(user_fnc)
          .sublist("function-exprtk")
          .set<int>("number of arguments", dim_ + 1)
          .set<std::string>("formula", bcs.formulas[i]);
        forms_copy[i] = user_fnc;
      }
    }
    bcfn.sublist("function-tabular")
      .set<Teuchos::Array<double>>("x values", bcs.times)
      .set<Teuchos::Array<double>>("y values", bcs.values)
      .set<Teuchos::Array<std::string>>("forms", forms_copy);
    flag = true;
  }

  return flag;
}

} // namespace AmanziInput
} // namespace Amanzi
