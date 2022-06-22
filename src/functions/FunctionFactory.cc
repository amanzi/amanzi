/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "errors.hh"
#include "HDF5Reader.hh"

#include "FunctionAdditive.hh"
#include "FunctionBilinear.hh"
#include "FunctionComposition.hh"
#include "FunctionConstant.hh"
#include "FunctionDistance.hh"
#include "FunctionSquareDistance.hh"
#include "FunctionFactory.hh"
#include "FunctionLinear.hh"
#include "FunctionMonomial.hh"
#include "FunctionMultiplicative.hh"
#include "FunctionPolynomial.hh"
#include "FunctionSeparable.hh"
#include "FunctionSmoothStep.hh"
#include "FunctionStandardMath.hh"
#include "FunctionStaticHead.hh"
#include "FunctionTabular.hh"

namespace Amanzi {

Function*
FunctionFactory::Create(Teuchos::ParameterList& list) const
{
  // Iterate through the parameters in the list.  There should be exactly
  // one, a sublist, whose name matches one of the known function types.
  // Anything else is a syntax error and we throw an exception.
  Function* f = 0;
  for (auto it = list.begin(); it != list.end(); ++it) {
    std::string function_type = list.name(it);
    if (list.isSublist(function_type)) { // process the function sublist
      if (f) { // error: already processed a function sublist
        Errors::Message m;
        m << "FunctionFactory: extraneous function sublist: "
          << function_type.c_str();
        Exceptions::amanzi_throw(m);
      }
      Teuchos::ParameterList& function_params = list.sublist(function_type);
      if (function_type == "function-constant")
        f = create_constant(function_params);
      else if (function_type == "function-tabular")
        f = create_tabular(function_params);
      else if (function_type == "function-polynomial")
        f = create_polynomial(function_params);
      else if (function_type == "function-monomial")
        f = create_monomial(function_params);
      else if (function_type == "function-smooth-step")
        f = create_smooth_step(function_params);
      else if (function_type == "function-linear")
        f = create_linear(function_params);
      else if (function_type == "function-separable")
        f = create_separable(function_params);
      else if (function_type == "function-additive")
        f = create_additive(function_params);
      else if (function_type == "function-multiplicative")
        f = create_multiplicative(function_params);
      else if (function_type == "function-composition")
        f = create_composition(function_params);
      else if (function_type == "function-static-head")
        f = create_static_head(function_params);
      else if (function_type == "function-standard-math")
        f = create_standard_math(function_params);
      else if (function_type == "function-bilinear")
        f = create_bilinear(function_params);
      else if (function_type == "function-distance")
        f = create_distance(function_params);
      else if (function_type == "function-squaredistance")
        f = create_squaredistance(function_params);
      else { // I don't recognize this function type
        if (f) delete f;
        Errors::Message m;
        m << "FunctionFactory: unknown function type: "
          << function_type.c_str();
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

Function*
FunctionFactory::create_constant(Teuchos::ParameterList& params) const
{
  Function* f;
  try {
    double value = params.get<double>("value");
    f = new FunctionConstant(value);
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-constant parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function*
FunctionFactory::create_tabular(Teuchos::ParameterList& params) const
{
  Function* f;

  if (params.isParameter("file")) {
    //    try {
    std::string filename = params.get<std::string>("file");
    HDF5Reader reader(filename);

    int xi = 0;
    std::string x = params.get<std::string>("x header");
    std::string xc = params.get<std::string>("x coordinate", "t");
    if (xc.compare(0, 1, "t") == 0)
      xi = 0;
    else if (xc.compare(0, 1, "x") == 0)
      xi = 1;
    else if (xc.compare(0, 1, "y") == 0)
      xi = 2;
    else if (xc.compare(0, 1, "z") == 0)
      xi = 3;
    std::string y = params.get<std::string>("y header");

    Kokkos::View<double*,Kokkos::HostSpace> vec_x;
    Kokkos::View<double*,Kokkos::HostSpace> vec_y;
    reader.ReadData(x, vec_x);
    reader.ReadData(y, vec_y);
    if (params.isParameter("forms")) {
      Teuchos::Array<std::string> form_strings(
        params.get<Teuchos::Array<std::string>>("forms"));
      Kokkos::View<FunctionTabular::Form*,Kokkos::HostSpace> form("form", form_strings.size());
      for (int i = 0; i < form_strings.size(); ++i) {
        if (form_strings[i] == "linear")
          form[i] = FunctionTabular::LINEAR;
        else if (form_strings[i] == "constant")
          form[i] = FunctionTabular::CONSTANT;
        else {
          Errors::Message m;
          m << "unknown form \"" << form_strings[i].c_str() << "\"";
          Exceptions::amanzi_throw(m);
        }
      }
      f = new FunctionTabular(vec_x, vec_y, xi, form);
    } else {
      f = new FunctionTabular(vec_x, vec_y, xi);
    }

    // }
    // catch (Teuchos::Exceptions::InvalidParameter& msg) {
    //   Errors::Message m;
    //   m << "FunctionFactory: function-tabular parameter error: " <<
    //   msg.what(); Exceptions::amanzi_throw(m);
    // }
    // catch (Errors::Message& msg) {
    //   Errors::Message m;
    //   m << "FunctionFactory: function-tabular parameter error: " <<
    //   msg.what(); Exceptions::amanzi_throw(m);
    // }
  } else {
    try {
      std::vector<double> x_vec(
        params.get<Teuchos::Array<double>>("x values").toVector());
      Kokkos::View<double*,Kokkos::HostSpace> x("x", x_vec.size());
      for (int i = 0; i < x.extent(0); ++i) x(i) = x_vec[i];

      std::string xc = params.get<std::string>("x coordinate", "t");
      int xi = 0;
      if (xc.compare(0, 1, "t") == 0)
        xi = 0;
      else if (xc.compare(0, 1, "x") == 0)
        xi = 1;
      else if (xc.compare(0, 1, "y") == 0)
        xi = 2;
      else if (xc.compare(0, 1, "z") == 0)
        xi = 3;

      std::vector<double> y_vec(
        params.get<Teuchos::Array<double>>("y values").toVector());
      Kokkos::View<double*,Kokkos::HostSpace> y("y", y_vec.size());
      for (int i = 0; i < y.extent(0); ++i) y(i) = y_vec[i];

      if (params.isParameter("forms")) {
        Teuchos::Array<std::string> form_strings(
          params.get<Teuchos::Array<std::string>>("forms"));
        int nforms = form_strings.size();
        Kokkos::View<FunctionTabular::Form*,Kokkos::HostSpace> form("form", nforms);

        bool flag_func(false);
        std::vector<Function*> func(nforms);

        for (int i = 0; i < nforms; ++i) {
          if (form_strings[i] == "linear")
            form[i] = FunctionTabular::LINEAR;
          else if (form_strings[i] == "constant")
            form[i] = FunctionTabular::CONSTANT;
          else {
            form[i] = FunctionTabular::FUNCTION;
            if (params.isSublist(form_strings[i])) {
              Teuchos::ParameterList& f1_params =
                params.sublist(form_strings[i]);

              Function* f1;
              FunctionFactory factory;
              f1 = factory.Create(f1_params);

              func[i] = f1;
              flag_func = true;
            } else {
              Errors::Message m;
              m << "unknown form \"" << form_strings[i].c_str() << "\"";
              Exceptions::amanzi_throw(m);
            }
          }
        }
        if (flag_func) {
          f = new FunctionTabular(x, y, xi, form, func);
        } else {
          f = new FunctionTabular(x, y, xi, form);
        }
      } else {
        f = new FunctionTabular(x, y, xi);
      }
    } catch (Teuchos::Exceptions::InvalidParameter& msg) {
      Errors::Message m;
      m << "FunctionFactory: function-tabular parameter error: " << msg.what();
      Exceptions::amanzi_throw(m);
    } catch (Errors::Message& msg) {
      Errors::Message m;
      m << "FunctionFactory: function-tabular parameter error: " << msg.what();
      Exceptions::amanzi_throw(m);
    }
  }
  return f;
}

Function*
FunctionFactory::create_smooth_step(Teuchos::ParameterList& params) const
{
  Function* f;
  try {
    double x0 = params.get<double>("x0");
    double x1 = params.get<double>("x1");
    double y0 = params.get<double>("y0");
    double y1 = params.get<double>("y1");
    f = new FunctionSmoothStep(x0, y0, x1, y1);
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-smooth-step parameter error: "
      << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-smooth-step parameter error: "
      << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function*
FunctionFactory::create_polynomial(Teuchos::ParameterList& params) const
{
  Function* f;
  try {
    std::vector<double> c_vec(
      params.get<Teuchos::Array<double>>("coefficients").toVector());
    std::vector<int> p_vec(
      params.get<Teuchos::Array<int>>("exponents").toVector());

    Kokkos::View<double*,Kokkos::HostSpace> c("c", c_vec.size());
    for (int i = 0; i < c.extent(0); ++i) c(i) = c_vec[i];
    Kokkos::View<int*,Kokkos::HostSpace> p("p", p_vec.size());
    for (int i = 0; i < p.extent(0); ++i) p(i) = p_vec[i];

    double x0 = params.get<double>("reference point", 0.0);
    f = new FunctionPolynomial(c, p, x0);
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-polynomial parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-polynomial parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function*
FunctionFactory::create_monomial(Teuchos::ParameterList& params) const
{
  Function* f;
  try {
    double c = params.get<double>("c");
    std::vector<double> x0_vec(
      params.get<Teuchos::Array<double>>("x0").toVector());
    std::vector<int> p_vec(
      params.get<Teuchos::Array<int>>("exponents").toVector());
    Kokkos::View<double*,Kokkos::HostSpace> x0("x0", x0_vec.size());
    for (int i = 0; i < x0.extent(0); ++i) x0(i) = x0_vec[i];
    Kokkos::View<int*,Kokkos::HostSpace> p("p", p_vec.size());
    for (int i = 0; i < p.extent(0); ++i) p(i) = p_vec[i];
    f = new FunctionMonomial(c, x0, p);
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-monomial parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-monomial parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function*
FunctionFactory::create_linear(Teuchos::ParameterList& params) const
{
  Function* f;
  try {
    double y0 = params.get<double>("y0");
    std::vector<double> grad_vec(
      params.get<Teuchos::Array<double>>("gradient").toVector());
    Teuchos::Array<double> zero(grad_vec.size(), 0.0);
    std::vector<double> x0_vec(
      params.get<Teuchos::Array<double>>("x0", zero).toVector());

    Kokkos::View<double*,Kokkos::HostSpace> grad("grad", grad_vec.size());
    for (int i = 0; i < grad.extent(0); ++i) grad(i) = grad_vec[i];
    Kokkos::View<double*,Kokkos::HostSpace> x0("x0", x0_vec.size());
    for (int i = 0; i < x0.extent(0); ++i) x0(i) = x0_vec[i];

    f = new FunctionLinear(y0, grad, x0);
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-linear parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-linear parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function*
FunctionFactory::create_separable(Teuchos::ParameterList& params) const
{
  Function* f;
  std::unique_ptr<Function> f1, f2;
  FunctionFactory factory;
  try {
    if (params.isSublist("function1")) {
      Teuchos::ParameterList& f1_params = params.sublist("function1");
      f1 = std::unique_ptr<Function>(factory.Create(f1_params));
    } else {
      Errors::Message m;
      m << "missing sublist function1";
      Exceptions::amanzi_throw(m);
    }
    if (params.isSublist("function2")) {
      Teuchos::ParameterList& f2_params = params.sublist("function2");
      f2 = std::unique_ptr<Function>(factory.Create(f2_params));
    } else {
      Errors::Message m;
      m << "missing sublist function2";
      Exceptions::amanzi_throw(m);
    }
    f = new FunctionSeparable(std::move(f1), std::move(f2));
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-separable parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function*
FunctionFactory::create_additive(Teuchos::ParameterList& params) const
{
  Function* f;
  std::unique_ptr<Function> f1, f2;
  FunctionFactory factory;
  try {
    if (params.isSublist("function1")) {
      Teuchos::ParameterList& f1_params = params.sublist("function1");
      f1 = std::unique_ptr<Function>(factory.Create(f1_params));
    } else {
      Errors::Message m;
      m << "missing sublist function1";
      Exceptions::amanzi_throw(m);
    }
    if (params.isSublist("function2")) {
      Teuchos::ParameterList& f2_params = params.sublist("function2");
      f2 = std::unique_ptr<Function>(factory.Create(f2_params));
    } else {
      Errors::Message m;
      m << "missing sublist function2";
      Exceptions::amanzi_throw(m);
    }
    f = new FunctionAdditive(std::move(f1), std::move(f2));
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-additive parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function*
FunctionFactory::create_multiplicative(Teuchos::ParameterList& params) const
{
  Function* f;
  std::unique_ptr<Function> f1, f2;
  FunctionFactory factory;
  try {
    if (params.isSublist("function1")) {
      Teuchos::ParameterList& f1_params = params.sublist("function1");
      f1 = std::unique_ptr<Function>(factory.Create(f1_params));
    } else {
      Errors::Message m;
      m << "missing sublist function1";
      Exceptions::amanzi_throw(m);
    }
    if (params.isSublist("function2")) {
      Teuchos::ParameterList& f2_params = params.sublist("function2");
      f2 = std::unique_ptr<Function>(factory.Create(f2_params));
    } else {
      Errors::Message m;
      m << "missing sublist function2";
      Exceptions::amanzi_throw(m);
    }
    f = new FunctionMultiplicative(std::move(f1), std::move(f2));
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-multiplicative parameter error: "
      << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function*
FunctionFactory::create_composition(Teuchos::ParameterList& params) const
{
  Function* f;
  std::unique_ptr<Function> f1, f2;
  FunctionFactory factory;
  try {
    if (params.isSublist("function1")) {
      Teuchos::ParameterList& f1_params = params.sublist("function1");
      f1 = std::unique_ptr<Function>(factory.Create(f1_params));
    } else {
      Errors::Message m;
      m << "missing sublist function1";
      Exceptions::amanzi_throw(m);
    }
    if (params.isSublist("function2")) {
      Teuchos::ParameterList& f2_params = params.sublist("function2");
      f2 = std::unique_ptr<Function>(factory.Create(f2_params));
    } else {
      Errors::Message m;
      m << "missing sublist function2";
      Exceptions::amanzi_throw(m);
    }
    f = new FunctionComposition(std::move(f1), std::move(f2));
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-composition parameter error: "
      << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function*
FunctionFactory::create_static_head(Teuchos::ParameterList& params) const
{
  Function* f = nullptr;
  FunctionFactory factory;
  try {
    double p0 = params.get<double>("p0");
    double density = params.get<double>("density");
    double gravity = params.get<double>("gravity");
    int dim = params.get<int>("space dimension");
    if (params.isSublist("water table elevation")) {
      Teuchos::ParameterList& sublist = params.sublist("water table elevation");
      std::unique_ptr<Function> water_table(factory.Create(sublist));
      f = new FunctionStaticHead(
        p0, density, gravity, std::move(water_table), dim);
    } else {
      Errors::Message m;
      m << "missing sublist \"water table elevation\"";
      Exceptions::amanzi_throw(m);
    }
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-static-head parameter error: "
      << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-static-head parameter error: "
      << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function*
FunctionFactory::create_standard_math(Teuchos::ParameterList& params) const
{
  Function* f;
  FunctionFactory factory;
  try {
    std::string op_string = params.get<std::string>("operator");
    char op[10];
    strcpy(op, op_string.c_str());
    double amplitude = params.get<double>("amplitude", 1.0);
    double param = params.get<double>("parameter", 1.0);
    double shift = params.get<double>("shift", 0.0);
    f = new FunctionStandardMath(op, amplitude, param, shift);
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-standard-math parameter error: "
      << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-standard-math parameter error: "
      << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function*
FunctionFactory::create_bilinear(Teuchos::ParameterList& params) const
{
  Function* f = nullptr;

  if (params.isParameter("file")) {
    try {
      std::string filename = params.get<std::string>("file");
      HDF5Reader reader(filename);

      int xi, yi(0); // input indices
      std::string x = params.get<std::string>("row header");
      std::string xdim = params.get<std::string>("row coordinate");
      if (xdim.compare(0, 1, "t") == 0)
        xi = 0;
      else if (xdim.compare(0, 1, "x") == 0)
        xi = 1;
      else if (xdim.compare(0, 1, "y") == 0)
        xi = 2;
      else if (xdim.compare(0, 1, "z") == 0)
        xi = 3;
      else {
        Errors::Message m;
        m << "FunctionFactory: function-bilinear parameter error: invalid "
             "\"row coordinate\" \""
          << xdim << "\" must be one of \"t,\" \"x,\" \"y,\" \"z.\"";
        Exceptions::amanzi_throw(m);
        xi = 0;
      }
      Kokkos::View<double*,Kokkos::HostSpace> vec_x;
      reader.ReadData(x, vec_x);

      std::string y = params.get<std::string>("column header");
      std::string ydim = params.get<std::string>("column coordinate");
      if (ydim.compare(0, 1, "t") == 0)
        yi = 0;
      else if (ydim.compare(0, 1, "x") == 0)
        yi = 1;
      else if (ydim.compare(0, 1, "y") == 0)
        yi = 2;
      else if (ydim.compare(0, 1, "z") == 0)
        yi = 3;
      else {
        Errors::Message m;
        m << "FunctionFactory: function-bilinear parameter error: invalid "
             "\"column coordinate\" \""
          << ydim << "\" must be one of \"t,\" \"x,\" \"y,\" \"z.\"";
        Exceptions::amanzi_throw(m);
        yi = 0;
      }
      Kokkos::View<double*,Kokkos::HostSpace> vec_y;
      reader.ReadData(y, vec_y);

      std::string v = params.get<std::string>("value header");
      Kokkos::View<double**,Kokkos::HostSpace> mat_v;
      reader.ReadMatData(v, mat_v);

      f = new FunctionBilinear(vec_x, vec_y, mat_v, xi, yi);
    } catch (Teuchos::Exceptions::InvalidParameter& msg) {
      Errors::Message m;
      m << "FunctionFactory: function-bilinear parameter error: " << msg.what();
      Exceptions::amanzi_throw(m);
    } catch (Errors::Message& msg) {
      Errors::Message m;
      m << "FunctionFactory: function-bilinear parameter error: " << msg.what();
      Exceptions::amanzi_throw(m);
    }
  } else {
    Errors::Message m;
    m << "missing parameter \"file\"";
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function*
FunctionFactory::create_distance(Teuchos::ParameterList& params) const
{
  Function* f;
  try {
    std::vector<double> x0_vec(
      params.get<Teuchos::Array<double>>("x0").toVector());
    std::vector<double> metric_vec(
      params.get<Teuchos::Array<double>>("metric").toVector());

    Kokkos::View<double*,Kokkos::HostSpace> x0("x0", x0_vec.size());
    for (int i = 0; i < x0.extent(0); ++i) x0(i) = x0_vec[i];
    Kokkos::View<double*,Kokkos::HostSpace> metric("metric", metric_vec.size());
    for (int i = 0; i < metric.extent(0); ++i) metric(i) = metric_vec[i];

    f = new FunctionDistance(x0, metric);
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-distance parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-distance parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

Function*
FunctionFactory::create_squaredistance(Teuchos::ParameterList& params) const
{
  Function* f;
  try {
    std::vector<double> x0_vec(
      params.get<Teuchos::Array<double>>("x0").toVector());
    std::vector<double> metric_vec(
      params.get<Teuchos::Array<double>>("metric").toVector());

    Kokkos::View<double*,Kokkos::HostSpace> x0("x0", x0_vec.size());
    for (int i = 0; i < x0.extent(0); ++i) x0(i) = x0_vec[i];
    Kokkos::View<double*,Kokkos::HostSpace> metric("metric", metric_vec.size());
    for (int i = 0; i < metric.extent(0); ++i) metric(i) = metric_vec[i];

    f = new FunctionSquareDistance(x0, metric);
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-squaredistance parameter error: "
      << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-squaredistance parameter error: "
      << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}
} // namespace Amanzi
