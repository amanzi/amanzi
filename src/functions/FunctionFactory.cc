/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "errors.hh"
#include "HDF5Reader.hh"
#include "Key.hh"

#include "FunctionAdditive.hh"
#include "FunctionBilinear.hh"
#include "FunctionComposition.hh"
#include "FunctionConstant.hh"
#include "FunctionDistance.hh"
#include "FunctionExprTK.hh"
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
#include "FunctionBilinearAndTime.hh"

namespace Amanzi {

std::unique_ptr<Function>
FunctionFactory::Create(const std::string& function_type,
                        Teuchos::ParameterList& function_params) const
{
  std::unique_ptr<Function> f;
  if (function_type == "constant")
    f = create_constant(function_params);
  else if (function_type == "tabular")
    f = create_tabular(function_params);
  else if (function_type == "polynomial")
    f = create_polynomial(function_params);
  else if (function_type == "monomial")
    f = create_monomial(function_params);
  else if (function_type == "smooth step")
    f = create_smooth_step(function_params);
  else if (function_type == "linear")
    f = create_linear(function_params);
  else if (function_type == "separable")
    f = create_separable(function_params);
  else if (function_type == "additive")
    f = create_additive(function_params);
  else if (function_type == "multiplicative")
    f = create_multiplicative(function_params);
  else if (function_type == "composition")
    f = create_composition(function_params);
  else if (function_type == "static head")
    f = create_static_head(function_params);
  else if (function_type == "standard math")
    f = create_standard_math(function_params);
  else if (function_type == "bilinear")
    f = create_bilinear(function_params);
  else if (function_type == "bilinear and time")
    f = create_bilinear_and_time(function_params);
  else if (function_type == "distance")
    f = create_distance(function_params);
  else if (function_type == "squaredistance")
    f = create_squaredistance(function_params);
  else if (function_type == "exprtk")
    f = create_exprtk(function_params);
  else { // I don't recognize this function type
    Errors::Message m;
    m << "FunctionFactory: unknown function type: " << function_type.c_str();
    Exceptions::amanzi_throw(m);
  }
  return f;
}


std::unique_ptr<Function>
FunctionFactory::Create(Teuchos::ParameterList& list) const
{
  std::unique_ptr<Function> f;

  if (list.isParameter("function type")) {
    // new standardization of accessing typed things
    auto function_type = list.get<std::string>("function type");
    f = Create(function_type, list);

  } else {
    // Old style, deprecate this see #181
    //
    // Iterate through the parameters in the list.  There should be exactly
    // one, a sublist, whose name matches one of the known function types.
    // Anything else is a syntax error and we throw an exception.
    for (auto it = list.begin(); it != list.end(); ++it) {
      std::string function_type = list.name(it);
      if (list.isSublist(function_type)) { // process the function sublist
        if (f.get()) {                     // error: already processed a function sublist
          Errors::Message m;
          m << "FunctionFactory: extraneous function sublist: " << function_type.c_str();
          Exceptions::amanzi_throw(m);
        }
      } else { // not the expected function sublist
        Errors::Message m;
        m << "FunctionFactory: unknown parameter: " << function_type.c_str();
        Exceptions::amanzi_throw(m);
      }
      // strip the function-
      if (!Keys::starts_with(function_type, "function-")) {
        Errors::Message m;
        m << "FunctionFactory: unknown function type: \"" << function_type << "\"";
        Exceptions::amanzi_throw(m);
      }

      Teuchos::ParameterList& function_params = list.sublist(function_type);
      function_type = function_type.substr(std::string("function-").size(), std::string::npos);
      function_type = Keys::replace_all(function_type, "-", " ");
      f = Create(function_type, function_params);
    }

    if (!f) { // no function sublist was found above
      Errors::Message m;
      m << "FunctionFactory: missing function sublist.";
      Exceptions::amanzi_throw(m);
    }
  }
  return f;
}


std::unique_ptr<Function>
FunctionFactory::create_constant(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  try {
    double value = params.get<double>("value");
    f = std::make_unique<FunctionConstant>(value);
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-constant parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

std::unique_ptr<Function>
FunctionFactory::create_tabular(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;

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

    Kokkos::View<double*, Kokkos::HostSpace> vec_x;
    Kokkos::View<double*, Kokkos::HostSpace> vec_y;
    reader.ReadData(x, vec_x);
    reader.ReadData(y, vec_y);
    if (params.isParameter("forms")) {
      Kokkos::View<Form_kind*, Kokkos::HostSpace> form("forms");

      if (params.isType<Teuchos::Array<std::string>>("forms")) {
        Teuchos::Array<std::string> form_strings(params.get<Teuchos::Array<std::string>>("forms"));
        Kokkos::resize(form, form_strings.size());
        for (int i = 0; i < form_strings.size(); ++i) {
          if (form_strings[i] == "linear") form[i] = Form_kind::LINEAR;
          else if (form_strings[i] == "constant") form[i] = Form_kind::CONSTANT;
          else {
            Errors::Message m;
            m << "unknown form \"" << form_strings[i] << "\"";
            Exceptions::amanzi_throw(m);
          }
        }

      } else if (params.isType<std::string>("forms")) {
        std::string form_string = params.get<std::string>("forms");
        Kokkos::resize(form, vec_x.size() - 1);

        if (form_string == "linear") {
          Kokkos::deep_copy(form, Form_kind::LINEAR);
        } else if (form_string == "constant") {
          Kokkos::deep_copy(form, Form_kind::CONSTANT);
        } else {
          Errors::Message m;
          m << "unknown form \"" << form_string << "\"";
          Exceptions::amanzi_throw(m);
        }
      }
      f = std::make_unique<FunctionTabular>(vec_x, vec_y, xi, form);
    } else {
      f = std::make_unique<FunctionTabular>(vec_x, vec_y, xi);
    }

  } else {
    try {
      std::vector<double> x_vec(params.get<Teuchos::Array<double>>("x values").toVector());
      Kokkos::View<double*, Kokkos::HostSpace> x("x", x_vec.size());
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

      std::vector<double> y_vec(params.get<Teuchos::Array<double>>("y values").toVector());
      Kokkos::View<double*, Kokkos::HostSpace> y("y", y_vec.size());
      for (int i = 0; i < y.extent(0); ++i) y(i) = y_vec[i];

      if (params.isParameter("forms")) {
        Teuchos::Array<std::string> form_strings(params.get<Teuchos::Array<std::string>>("forms"));
        int nforms = form_strings.size();
        Kokkos::View<Form_kind*, Kokkos::HostSpace> form("form", nforms);

        bool flag_func(false);
        std::vector<std::unique_ptr<Function>> func(nforms);

        for (int i = 0; i < nforms; ++i) {
          if (form_strings[i] == "linear")
            form[i] = Form_kind::LINEAR;
          else if (form_strings[i] == "constant")
            form[i] = Form_kind::CONSTANT;
          else {
            form[i] = Form_kind::FUNCTION;
            if (params.isSublist(form_strings[i])) {
              Teuchos::ParameterList& f1_params = params.sublist(form_strings[i]);

              std::unique_ptr<Function> f1;
              FunctionFactory factory;
              f1 = factory.Create(f1_params);

              func[i] = std::move(f1);
              flag_func = true;
            } else {
              Errors::Message m;
              m << "unknown form \"" << form_strings[i].c_str() << "\"";
              Exceptions::amanzi_throw(m);
            }
          }
        }
        if (flag_func) {
          f = std::make_unique<FunctionTabular>(x, y, xi, form, std::move(func));
        } else {
          f = std::make_unique<FunctionTabular>(x, y, xi, form);
        }
      } else {
        f = std::make_unique<FunctionTabular>(x, y, xi);
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

std::unique_ptr<Function>
FunctionFactory::create_smooth_step(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  try {
    double x0 = params.get<double>("x0");
    double x1 = params.get<double>("x1");
    double y0 = params.get<double>("y0");
    double y1 = params.get<double>("y1");
    f = std::make_unique<FunctionSmoothStep>(x0, y0, x1, y1);
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-smooth-step parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-smooth-step parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

std::unique_ptr<Function>
FunctionFactory::create_polynomial(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  try {
    std::vector<double> c_vec(params.get<Teuchos::Array<double>>("coefficients").toVector());
    std::vector<int> p_vec(params.get<Teuchos::Array<int>>("exponents").toVector());

    Kokkos::View<double*, Kokkos::HostSpace> c("c", c_vec.size());
    for (int i = 0; i < c.extent(0); ++i) c(i) = c_vec[i];
    Kokkos::View<int*, Kokkos::HostSpace> p("p", p_vec.size());
    for (int i = 0; i < p.extent(0); ++i) p(i) = p_vec[i];

    double x0 = params.get<double>("reference point", 0.0);
    f = std::make_unique<FunctionPolynomial>(c, p, x0);
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

std::unique_ptr<Function>
FunctionFactory::create_monomial(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  try {
    double c = params.get<double>("c");
    std::vector<double> x0_vec(params.get<Teuchos::Array<double>>("x0").toVector());
    std::vector<int> p_vec(params.get<Teuchos::Array<int>>("exponents").toVector());
    Kokkos::View<double*, Kokkos::HostSpace> x0("x0", x0_vec.size());
    for (int i = 0; i < x0.extent(0); ++i) x0(i) = x0_vec[i];
    Kokkos::View<int*, Kokkos::HostSpace> p("p", p_vec.size());
    for (int i = 0; i < p.extent(0); ++i) p(i) = p_vec[i];
    f = std::make_unique<FunctionMonomial>(c, x0, p);
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

std::unique_ptr<Function>
FunctionFactory::create_linear(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  try {
    double y0 = params.get<double>("y0");
    std::vector<double> grad_vec(params.get<Teuchos::Array<double>>("gradient").toVector());
    Teuchos::Array<double> zero(grad_vec.size(), 0.0);
    std::vector<double> x0_vec(params.get<Teuchos::Array<double>>("x0", zero).toVector());

    Kokkos::View<double*, Kokkos::HostSpace> grad("grad", grad_vec.size());
    for (int i = 0; i < grad.extent(0); ++i) grad(i) = grad_vec[i];
    Kokkos::View<double*, Kokkos::HostSpace> x0("x0", x0_vec.size());
    for (int i = 0; i < x0.extent(0); ++i) x0(i) = x0_vec[i];

    f = std::make_unique<FunctionLinear>(y0, grad, x0);
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

std::unique_ptr<Function>
FunctionFactory::create_separable(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  std::unique_ptr<Function> f1, f2;
  FunctionFactory factory;
  try {
    if (params.isSublist("function1")) {
      Teuchos::ParameterList& f1_params = params.sublist("function1");
      f1 = factory.Create(f1_params);
    } else {
      Errors::Message m;
      m << "missing sublist function1";
      Exceptions::amanzi_throw(m);
    }
    if (params.isSublist("function2")) {
      Teuchos::ParameterList& f2_params = params.sublist("function2");
      f2 = factory.Create(f2_params);
    } else {
      Errors::Message m;
      m << "missing sublist function2";
      Exceptions::amanzi_throw(m);
    }
    f = std::make_unique<FunctionSeparable>(std::move(f1), std::move(f2));
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-separable parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

std::unique_ptr<Function>
FunctionFactory::create_additive(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  std::unique_ptr<Function> f1, f2;
  FunctionFactory factory;
  try {
    if (params.isSublist("function1")) {
      Teuchos::ParameterList& f1_params = params.sublist("function1");
      f1 = factory.Create(f1_params);
    } else {
      Errors::Message m;
      m << "missing sublist function1";
      Exceptions::amanzi_throw(m);
    }
    if (params.isSublist("function2")) {
      Teuchos::ParameterList& f2_params = params.sublist("function2");
      f2 = factory.Create(f2_params);
    } else {
      Errors::Message m;
      m << "missing sublist function2";
      Exceptions::amanzi_throw(m);
    }
    f = std::make_unique<FunctionAdditive>(std::move(f1), std::move(f2));
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-additive parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

std::unique_ptr<Function>
FunctionFactory::create_multiplicative(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  std::unique_ptr<Function> f1, f2;
  FunctionFactory factory;
  try {
    if (params.isSublist("function1")) {
      Teuchos::ParameterList& f1_params = params.sublist("function1");
      f1 = factory.Create(f1_params);
    } else {
      Errors::Message m;
      m << "missing sublist function1";
      Exceptions::amanzi_throw(m);
    }
    if (params.isSublist("function2")) {
      Teuchos::ParameterList& f2_params = params.sublist("function2");
      f2 = factory.Create(f2_params);
    } else {
      Errors::Message m;
      m << "missing sublist function2";
      Exceptions::amanzi_throw(m);
    }
    f = std::make_unique<FunctionMultiplicative>(std::move(f1), std::move(f2));
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-multiplicative parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

std::unique_ptr<Function>
FunctionFactory::create_composition(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  std::unique_ptr<Function> f1, f2;
  FunctionFactory factory;
  try {
    if (params.isSublist("function1")) {
      Teuchos::ParameterList& f1_params = params.sublist("function1");
      f1 = factory.Create(f1_params);
    } else {
      Errors::Message m;
      m << "missing sublist function1";
      Exceptions::amanzi_throw(m);
    }
    if (params.isSublist("function2")) {
      Teuchos::ParameterList& f2_params = params.sublist("function2");
      f2 = factory.Create(f2_params);
    } else {
      Errors::Message m;
      m << "missing sublist function2";
      Exceptions::amanzi_throw(m);
    }
    f = std::make_unique<FunctionComposition>(std::move(f1), std::move(f2));
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-composition parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

std::unique_ptr<Function>
FunctionFactory::create_static_head(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  FunctionFactory factory;
  try {
    double p0 = params.get<double>("p0");
    double density = params.get<double>("density");
    double gravity = params.get<double>("gravity");
    int dim = params.get<int>("space dimension");
    if (params.isSublist("water table elevation")) {
      Teuchos::ParameterList& sublist = params.sublist("water table elevation");
      auto water_table = factory.Create(sublist);
      f = std::make_unique<FunctionStaticHead>(p0, density, gravity, std::move(water_table), dim);
    } else {
      Errors::Message m;
      m << "missing sublist \"water table elevation\"";
      Exceptions::amanzi_throw(m);
    }
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-static-head parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-static-head parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

std::unique_ptr<Function>
FunctionFactory::create_standard_math(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  FunctionFactory factory;
  try {
    std::string op_string = params.get<std::string>("operator");
    char op[10];
    strcpy(op, op_string.c_str());
    double amplitude = params.get<double>("amplitude", 1.0);
    double param = params.get<double>("parameter", 1.0);
    double shift = params.get<double>("shift", 0.0);
    f = std::make_unique<FunctionStandardMath>(op, amplitude, param, shift);
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-standard-math parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-standard-math parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

std::unique_ptr<Function>
FunctionFactory::create_bilinear(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;

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
      Kokkos::View<double*, Kokkos::HostSpace> vec_x;
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
      Kokkos::View<double*, Kokkos::HostSpace> vec_y;
      reader.ReadData(y, vec_y);

      std::string v = params.get<std::string>("value header");
      Kokkos::View<double**, Kokkos::HostSpace> mat_v;
      reader.ReadMatData(v, mat_v);
      f = std::make_unique<FunctionBilinear>(vec_x, vec_y, mat_v, xi, yi);
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

std::unique_ptr<Function>
FunctionFactory::create_distance(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  try {
    std::vector<double> x0_vec(params.get<Teuchos::Array<double>>("x0").toVector());
    std::vector<double> metric_vec(params.get<Teuchos::Array<double>>("metric").toVector());

    Kokkos::View<double*, Kokkos::HostSpace> x0("x0", x0_vec.size());
    for (int i = 0; i < x0.extent(0); ++i) x0(i) = x0_vec[i];
    Kokkos::View<double*, Kokkos::HostSpace> metric("metric", metric_vec.size());
    for (int i = 0; i < metric.extent(0); ++i) metric(i) = metric_vec[i];

    f = std::make_unique<FunctionDistance>(x0, metric);
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

std::unique_ptr<Function>
FunctionFactory::create_squaredistance(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  try {
    std::vector<double> x0_vec(params.get<Teuchos::Array<double>>("x0").toVector());
    std::vector<double> metric_vec(params.get<Teuchos::Array<double>>("metric").toVector());

    Kokkos::View<double*, Kokkos::HostSpace> x0("x0", x0_vec.size());
    for (int i = 0; i < x0.extent(0); ++i) x0(i) = x0_vec[i];
    Kokkos::View<double*, Kokkos::HostSpace> metric("metric", metric_vec.size());
    for (int i = 0; i < metric.extent(0); ++i) metric(i) = metric_vec[i];

    f = std::make_unique<FunctionSquareDistance>(x0, metric);
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-squaredistance parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-squaredistance parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}


std::unique_ptr<Function>
FunctionFactory::create_bilinear_and_time(Teuchos::ParameterList& params) const
{
  std::string filename = params.get<std::string>("file");
  std::string row_header = params.get<std::string>("row header");
  std::string row_coordinate = params.get<std::string>("row coordinate", "x");
  if (row_coordinate != "x" && row_coordinate != "y" && row_coordinate != "z") {
    Errors::Message m("FunctionFactory: function-bilinear_and_time does not do bilinear in time.");
    Exceptions::amanzi_throw(m);
  }
  std::string col_header = params.get<std::string>("column header");
  std::string col_coordinate = params.get<std::string>("column coordinate", "y");
  if (col_coordinate != "x" && col_coordinate != "y" && col_coordinate != "z") {
    Errors::Message m("FunctionFactory: function-bilinear_and_time does not do bilinear in time.");
    Exceptions::amanzi_throw(m);
  }
  std::string time_header = params.get<std::string>("time header");
  std::string value_header = params.get<std::string>("value header");

  std::string form_string = params.get<std::string>("forms", "linear");
  Form_kind form;
  if (form_string == "linear") {
    form = Form_kind::LINEAR;
  } else if (form_string == "constant") {
    form = Form_kind::CONSTANT;
  } else {
    Errors::Message m;
    m << "FunctionFactory: function-bilinear-and-time provided invalid value for \"form\" of \""
      << form_string << "\" -- valid are \"linear\" and \"constant\"";
    Exceptions::amanzi_throw(m);
  }
  return std::make_unique<FunctionBilinearAndTime>(filename,
                                                   time_header,
                                                   row_header,
                                                   row_coordinate,
                                                   col_header,
                                                   col_coordinate,
                                                   value_header,
                                                   form);
}


std::unique_ptr<Function>
FunctionFactory::create_exprtk(Teuchos::ParameterList& params) const
{
  std::unique_ptr<Function> f;
  try {
    int n = params.get<int>("number of arguments");
    std::string formula = params.get<std::string>("formula");
    f = std::make_unique<FunctionExprTK>(n, formula);
  } catch (Teuchos::Exceptions::InvalidParameter& msg) {
    Errors::Message m;
    m << "FunctionFactory: function-exprtk parameter error: " << msg.what();
    Exceptions::amanzi_throw(m);
  }
  return f;
}

} // namespace Amanzi
