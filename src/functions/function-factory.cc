#include "Teuchos_ParameterList.hpp"

#include "function-factory.hh"
#include "constant-function.hh"
#include "tabular-function.hh"
#include "smooth-step-function.hh"
#include "polynomial-function.hh"
#include "linear-function.hh"
#include "separable-function.hh"
#include "static-head-function.hh"
#include "errors.hh"

namespace Amanzi {

Function* FunctionFactory::Create(Teuchos::ParameterList &list) const
{
  // Iterate through the parameters in the list.  There should be exactly
  // one, a sublist, whose name matches one of the known function types.
  // Anything else is a syntax error and we throw an exception.
  Function *f = 0;
  for (Teuchos::ParameterList::ConstIterator it = list.begin(); it != list.end(); ++it) {
    std::string function_type = list.name(it);
    if (list.isSublist(function_type)) { // process the function sublist
      if (f) { // error: already processed a function sublist
        Errors::Message m;
        m << "FunctionFactory: extraneous function sublist: " << function_type.c_str();
        Exceptions::amanzi_throw(m);
      }
      Teuchos::ParameterList &function_params = list.sublist(function_type);
      if (function_type == "function-constant")
        f = create_constant(function_params);
      else if (function_type == "function-tabular")
        f = create_tabular(function_params);
      else if (function_type == "function-polynomial")
        f = create_polynomial(function_params);
      else if (function_type == "function-smooth-step")
        f = create_smooth_step(function_params);
      else if (function_type == "function-linear")
        f = create_linear(function_params);
      else if (function_type == "function-separable")
        f = create_separable(function_params);
      else if (function_type == "function-static-head")
        f = create_static_head(function_params);
      else {  // I don't recognize this function type
        if (f) delete f;
        Errors::Message m;
        m << "FunctionFactory: unknown function type: " << function_type.c_str();
        Exceptions::amanzi_throw(m);
      }
    } else { // not the expected function sublist
      if (f) delete f;
      Errors::Message m;
      m << "FunctionFactory: unknown parameter: " << function_type.c_str();
      Exceptions::amanzi_throw(m);
    }
  }
 
  if (!f) { // no function sublist was found above
    Errors::Message m;
    m << "FunctionFactory: missing function sublist.";
    Exceptions::amanzi_throw(m);
  }

  return f;
}

Function* FunctionFactory::create_constant(Teuchos::ParameterList &params) const
{
  Function *f;
  try {
    double value = params.get<double>("value");
    f = new ConstantFunction(value);
  } catch (Teuchos::Exceptions::InvalidParameter &msg) {
    Errors::Message m;
    m << "FunctionFactory: function-constant parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function* FunctionFactory::create_tabular(Teuchos::ParameterList &params) const
{
  Function *f;
  try {
    std::vector<double> x(params.get<Teuchos::Array<double> >("x values").toVector());
    std::vector<double> y(params.get<Teuchos::Array<double> >("y values").toVector());
    if (params.isParameter("forms")) {
      Teuchos::Array<std::string> form_strings(params.get<Teuchos::Array<std::string> >("forms"));
      std::vector<TabularFunction::Form> form(form_strings.size());
      for (int i = 0; i < form_strings.size(); ++i) {
        if (form_strings[i] == "linear")
          form[i] = TabularFunction::LINEAR;
        else if (form_strings[i] == "constant")
          form[i] = TabularFunction::CONSTANT;
        else {
          Errors::Message m;
          m << "unknown form \"" << form_strings[i].c_str() << "\"";
          Exceptions::amanzi_throw(m);
        }
      }
      f = new TabularFunction(x, y, form);
    } else {
      f = new TabularFunction(x, y);
    }
  }
  catch (Teuchos::Exceptions::InvalidParameter &msg) {
    Errors::Message m;
    m << "FunctionFactory: function-tabular parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  catch (Errors::Message &msg) {
    Errors::Message m;
    m << "FunctionFactory: function-tabular parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function* FunctionFactory::create_smooth_step(Teuchos::ParameterList &params) const
{
  Function *f;
  try {
    double x0 = params.get<double>("x0");
    double x1 = params.get<double>("x1");
    double y0 = params.get<double>("y0");
    double y1 = params.get<double>("y1");
    f = new SmoothStepFunction(x0, y0, x1, y1);
  } catch (Teuchos::Exceptions::InvalidParameter &msg) {
    Errors::Message m;
    m << "FunctionFactory: function-smooth-step parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message &msg) {
    Errors::Message m;
    m << "FunctionFactory: function-smooth-step parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function* FunctionFactory::create_polynomial(Teuchos::ParameterList &params) const
{
  Function *f;
  try {
    std::vector<double> c(params.get<Teuchos::Array<double> >("coefficients").toVector());
    std::vector<int> p(params.get<Teuchos::Array<int> >("exponents").toVector());
    double x0 = params.get<double>("reference point", 0.0);
    f = new PolynomialFunction(c, p, x0);
  } catch (Teuchos::Exceptions::InvalidParameter &msg) {
    Errors::Message m;
    m << "FunctionFactory: function-polynomial parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  catch (Errors::Message &msg) {
    Errors::Message m;
    m << "FunctionFactory: function-polynomial parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function* FunctionFactory::create_linear(Teuchos::ParameterList &params) const
{
  Function *f;
  try {
    double y0 = params.get<double>("y0");
    std::vector<double> grad(params.get<Teuchos::Array<double> >("gradient").toVector());
    Teuchos::Array<double> zero(grad.size(),0.0);
    std::vector<double> x0(params.get<Teuchos::Array<double> >("x0", zero).toVector());
    f = new LinearFunction(y0, grad, x0);
  } catch (Teuchos::Exceptions::InvalidParameter &msg) {
    Errors::Message m;
    m << "FunctionFactory: function-linear parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  catch (Errors::Message &msg) {
    Errors::Message m;
    m << "FunctionFactory: function-linear parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function* FunctionFactory::create_separable(Teuchos::ParameterList &params) const
{
  Function *f;
  std::auto_ptr<Function> f1, f2;
  FunctionFactory factory;
  try {
    if (params.isSublist("function1")) {
      Teuchos::ParameterList &f1_params = params.sublist("function1");
      f1 = std::auto_ptr<Function>(factory.Create(f1_params));
    } else {
      Errors::Message m;
      m << "missing sublist function1";
      Exceptions::amanzi_throw(m);
    }
    if (params.isSublist("function2")) {
      Teuchos::ParameterList &f2_params = params.sublist("function2");
      f2 = std::auto_ptr<Function>(factory.Create(f2_params));
    } else {
      Errors::Message m;
      m << "missing sublist function2";
      Exceptions::amanzi_throw(m);
    }
    f = new SeparableFunction(f1, f2);
  } catch (Errors::Message &msg) {
    Errors::Message m;
    m << "FunctionFactory: function-separable parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function* FunctionFactory::create_static_head(Teuchos::ParameterList &params) const
{
  Function *f;
  FunctionFactory factory;
  try {
    double p0 = params.get<double>("p0");
    double density = params.get<double>("density");
    double gravity = params.get<double>("gravity");
    if (params.isSublist("water table elevation")) {
      Teuchos::ParameterList &sublist = params.sublist("water table elevation");
      std::auto_ptr<Function> water_table(factory.Create(sublist));
      f = new StaticHeadFunction(p0, density, gravity, water_table);
    } else {
      Errors::Message m;
      m << "missing sublist \"water table elevation\"";
      Exceptions::amanzi_throw(m);
    }
  } catch (Teuchos::Exceptions::InvalidParameter &msg) {
    Errors::Message m;
    m << "FunctionFactory: function-static-head parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message &msg) {
    Errors::Message m;
    m << "FunctionFactory: function-static-head parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

} // namespace Amanzi
