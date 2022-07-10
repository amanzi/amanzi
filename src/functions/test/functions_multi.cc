/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include <map>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AmanziComm.hh"
#include "MultiFunction.hh"
#include "errors.hh"

#include "FunctionAdditive.hh"
#include "FunctionMultiplicative.hh"
#include "FunctionComposition.hh"
#include "FunctionSeparable.hh"

#include "FunctionBilinear.hh"
#include "FunctionConstant.hh"
#include "FunctionDistance.hh"
#include "FunctionLinear.hh"
#include "FunctionMonomial.hh"
#include "FunctionPolynomial.hh"
#include "FunctionSmoothStep.hh"
#include "FunctionSquareDistance.hh"
#include "FunctionStandardMath.hh"
#include "FunctionStaticHead.hh"
#include "FunctionTabular.hh"

using namespace Amanzi;

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize(argc, argv);
  auto status = UnitTest::RunAllTests();
  Kokkos::finalize();
  return status;
}


struct test {};

// Test all the functions with both host and device
TEST_FIXTURE(test, values1)
{
  srand(20);
  constexpr size_t dims = 3;
  constexpr size_t nvalues = 20;
  size_t nfunctions = 0;
  Kokkos::DualView<double**> xyz_nolayout("xyz_nolayout", dims, nvalues);
  for (int i = 0; i < xyz_nolayout.extent(0); ++i) {
    for (int j = 0; j < xyz_nolayout.extent(1); ++j) {
      // Values between [0-10]
      xyz_nolayout.view_host()(i, j) =
        (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * 10.;
    }
  }
  Kokkos::deep_copy(xyz_nolayout.view_device(),xyz_nolayout.view_host()); 

  Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace> xyz_layout_left_h(
    "xyz_layout_left_h", dims, nvalues);
  for (int i = 0; i < xyz_layout_left_h.extent(0); ++i) {
    for (int j = 0; j < xyz_layout_left_h.extent(1); ++j) {
      // Values between [0-10]
      xyz_layout_left_h(i,j) = xyz_nolayout.view_host()(i, j);
    }
  }

  // vector of functions
  std::vector<Teuchos::RCP<const Function>> functions;

  // Function Constant --------------------------------------------
  functions.push_back(Teuchos::rcp(new FunctionConstant(1.5)));
 
  // Function Bilinear --------------------------------------------
  Kokkos::View<double*,Kokkos::HostSpace> x_b("x_b", 4);
  Kokkos::View<double*,Kokkos::HostSpace> y_b("y_b", 4);
  Kokkos::View<double**,Kokkos::HostSpace> v_b("v_b", 4, 4);
  x_b(0) = 0.0;
  x_b(1) = 2;
  x_b(2) = 4;
  x_b(3) = 8;
  y_b(0) = 2;
  y_b(1) = 4;
  y_b(2) = 6;
  y_b(3) = 10;
  v_b(0, 0) = 1.0;
  v_b(0, 1) = 3.5;
  v_b(0, 2) = 7.0;
  v_b(0, 3) = 9.5;
  v_b(1, 0) = 1.0;
  v_b(1, 1) = 3.5;
  v_b(1, 2) = 7.0;
  v_b(1, 3) = 9.5;
  v_b(2, 0) = 1.0;
  v_b(2, 1) = 3.5;
  v_b(2, 2) = 7.0;
  v_b(2, 3) = 9.5;
  v_b(3, 0) = 1.0;
  v_b(3, 1) = 3.5;
  v_b(3, 2) = 7.0;
  v_b(3, 3) = 9.5;
  functions.push_back(Teuchos::rcp(new FunctionBilinear(x_b, y_b, v_b, 0, 0)));
 
  // Function Distance -------------------------------------------
  Kokkos::View<double*,Kokkos::HostSpace> x0("x_a", dims);
  Kokkos::View<double*,Kokkos::HostSpace> metric("metric", dims);
  metric(0) = 1;
  metric(1) = 1;
  metric(2) = 1;
  x0(0) = xyz_nolayout.view_host()(0, 0);
  x0(1) = xyz_nolayout.view_host()(1, 0);
  x0(2) = xyz_nolayout.view_host()(2, 0);
  functions.push_back(Teuchos::rcp(new FunctionDistance(x0, metric)));

 
  // Function Linear ---------------------------------------------
  Kokkos::View<double*,Kokkos::HostSpace> grad("grad", dims);
  grad(0) = 1;
  grad(1) = 2;
  grad(2) = 3;
  functions.push_back(Teuchos::rcp(new FunctionLinear(0, grad)));

  // Function Monomial -------------------------------------------
  Kokkos::View<int*,Kokkos::HostSpace> p("p", dims);
  p(0) = 1;
  p(1) = 2;
  p(2) = 3;
  functions.push_back(Teuchos::rcp(new FunctionMonomial(1.25, x0, p)));

  // Function Polynomial------------------------------------------
  Kokkos::View<double*,Kokkos::HostSpace> c("c", dims);
  c(0) = 1;
  c(1) = 2;
  c(2) = 3;
  functions.push_back(Teuchos::rcp(new FunctionPolynomial(c, p)));

  // Function Additive  ------------------------------------------
  std::unique_ptr<Function> f1(new FunctionPolynomial(c, p));
  std::unique_ptr<Function> f2(new FunctionMonomial(1., x0, p));
  functions.push_back(
    Teuchos::rcp(new FunctionAdditive(std::move(f1), std::move(f2))));

  // Function Multiplicative -------------------------------------
  std::unique_ptr<Function> f3(new FunctionConstant(0.01));
  std::unique_ptr<Function> f4(new FunctionMonomial(1., x0, p));
  functions.push_back(
    Teuchos::rcp(new FunctionMultiplicative(std::move(f3), std::move(f4))));

  // Function SmoothStep -----------------------------------------
  functions.push_back(Teuchos::rcp(new FunctionSmoothStep(1, 2, 3, 4)));

  // Function SquareDistance -------------------------------------
  functions.push_back(Teuchos::rcp(new FunctionSquareDistance(x0, metric)));

  // Function standard maths -------------------------------------
  functions.push_back(
    Teuchos::rcp(new FunctionStandardMath("cos", 1., 1., 0.)));
  functions.push_back(
    Teuchos::rcp(new FunctionStandardMath("log", 1., 1., 0.)));
  functions.push_back(
    Teuchos::rcp(new FunctionStandardMath("pow", 1., 1., 0.)));
  functions.push_back(
    Teuchos::rcp(new FunctionStandardMath("ceil", 1., 1., 0.)));

  // Function Composition ----------------------------------------
  std::unique_ptr<Function> f5(new FunctionConstant(0.01));
  std::unique_ptr<Function> f6(new FunctionMonomial(1., x0, p));
  functions.push_back(
    Teuchos::rcp(new FunctionComposition(std::move(f5), std::move(f6))));

  // Function Separable ------------------------------------------
  std::unique_ptr<Function> f7(new FunctionConstant(0.01));
  Kokkos::View<double*,Kokkos::HostSpace> x0_l("x0_l", dims - 1);
  Kokkos::View<int*,Kokkos::HostSpace> p_l("p_l", dims - 1);
  p_l(0) = 1;
  p_l(1) = 2;
  x0_l(0) = xyz_nolayout.view_host()(0, 0);
  x0_l(1) = xyz_nolayout.view_host()(1, 0);
  std::unique_ptr<Function> f8(new FunctionMonomial(1., x0_l, p_l));
  functions.push_back(
    Teuchos::rcp(new FunctionSeparable(std::move(f7), std::move(f8))));

  // Function StaticHead -----------------------------------------
  std::unique_ptr<Function> f9(new FunctionMonomial(1., x0, p));
  functions.push_back(
    Teuchos::rcp(new FunctionStaticHead(1.0, 1.0, 1.0, std::move(f9), 2)));

  // Function Tabular --------------------------------------------
  Kokkos::View<double*,Kokkos::HostSpace> x0_i("x0_i", dims);
  x0_i(0) = 1.0;
  x0_i(1) = 2.0;
  x0_i(2) = 3.0;
  functions.push_back(Teuchos::rcp(new FunctionTabular(x0_i, x0_i, 0)));

  nfunctions = functions.size();

  // Create the multifunction
  MultiFunction mf = MultiFunction(functions);

  Kokkos::View<double**,Kokkos::HostSpace> result("result", nvalues, nfunctions);
  Kokkos::View<double**, Kokkos::LayoutLeft> result_gpu(
    "result", nvalues, nfunctions);

  mf.apply(xyz_nolayout.view_device(), result_gpu);

  for (int i = 0; i < nvalues; ++i) {
    Kokkos::View<double*,Kokkos::HostSpace> t = Kokkos::subview(xyz_layout_left_h, Kokkos::ALL, i);
    Kokkos::View<double*,Kokkos::HostSpace> tmp_res = mf(t);
    for (int j = 0; j < tmp_res.extent(0); ++j) { result(i, j) = tmp_res(j); }
  }
  // Copy GPU computation back 
  Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace> result_gpu_host("result gpu host",nvalues,nfunctions);
  Kokkos::deep_copy(result_gpu_host,result_gpu);  

  for (int j = 0; j < nfunctions; ++j) {
    for (int i = 0; i < nvalues; ++i) {
      if (fabs(result_gpu_host(i, j) - result(i, j)) >= 1.0e-10) {
        std::cerr << "Function: " << j << " ";
        std::cerr << result(i, j) << " == ";
        std::cerr << result_gpu_host(i, j) << " ; ";
      }
      CHECK_CLOSE(result(i, j), result_gpu_host(i, j), 1.0e-10);
    }
  }
 
}
