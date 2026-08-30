/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Energy PK
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CommonDefs.hh"
#include "CompositeVector.hh"
#include "evaluators_reg.hh"
#include "IAPWS95.hh"
#include "IAPWS95_RaggedSplinePH.hh"
#include "MeshFactory.hh"

TEST(ENTROPY_IAPWS95_SPLINE_DERIVATIVES)
{
  using namespace Amanzi;
  std::cout << "Test: Entropy spline properties" << std::endl;

  Teuchos::ParameterList plist;
  AmanziEOS::IAPWS95 eos95(plist);
  AmanziEOS::IAPWS95_RaggedSplinePH eos_spline(plist);

  eos_spline.InitializeSharedData();

  int n = 100;
  double eps(2e-6);
  double dp(0.1), dh(2.2), hl, hr, der;

  const auto& opt = eos_spline.GetOptions();
  for (double p = opt.p_min + dp; p < opt.p_max - dp; p += dp) {
    for (double h = opt.h_min + dh; h < opt.h_max - dh; h += dh) {
      const AmanziEOS::SaturationState& sat = eos95.SaturationLineP(p);
std::cout << p << " " << h << std::endl;
      if (eos_spline.IsPhysical(p, h, sat)) {
std::cout << "in: " << p << " " << h << std::endl;
        hl = h * (1.0 - eps);
        hr = h * (1.0 + eps);

        auto [prop1, liquid1, vapor1] = eos_spline.ThermodynamicsPH(p, h);
        auto [prop3, liquid3, vapor3] = eos_spline.ThermodynamicsPH(p, hr);
        auto [prop4, liquid4, vapor4] = eos_spline.ThermodynamicsPH(p, hl);

        // Second-order entropy derivative with respect to h 
        der = (prop4.s - 2 * prop1.s + prop3.s) / (eps * eps * h * h);
        const auto& exact = eos_spline.EntropyDerivativesPH(p, h);
        CHECK_CLOSE(der, exact[5], eps * std::fabs(exact[5]));
      }
    }
  }
}


/*
TEST(HELMHOLTZ_IAPWS95_SPLINE_EFFICIENCY)
{
  using namespace Amanzi;
  std::cout << "Test: Helmholtz energy spline efficiency" << std::endl;

  Teuchos::ParameterList plist;
  AmanziEOS::IAPWS95 eos95(plist);
  AmanziEOS::IAPWS95_RaggedSplineRhoT spline(plist);

  spline.InitializeSharedData();

  const int n = 100;
  double R(spline.R), dp(0.1), dT(2.2), p, T;

  std::ofstream out("eos_iapws95_spline.dat");
  for (int i = -n; i < n; ++i) {
    for (int j = -n; j < n; ++j) {
      p = eos95.PC + i * dp;
      T = eos95.TC + j * dT;

      // auto [prop, liquid, vapor] = eos95.ThermodynamicsPT(p, T);
      auto [prop, liquid, vapor] = spline.ThermodynamicsPT(p, T);
      // out << T << " " << p << " " << prop.helmholtz / R / T << "\n";
    }
  }
  out.close();

  auto calls = spline.residual_calls;
  auto itrs = spline.brent_root_itrs;
  double mean_calls = (double)calls / (4 * n * n);
  double mean_itrs = (double)itrs / (4 * n * n);
  std::cout << "Residual evaluations: " << calls << " or " << mean_calls << " times per sample\n";
  std::cout << "Brent's root find itrs: " << itrs << " or " << mean_itrs << " iterations per sample\n";
  std::cout << "Brent's bracket find itrs: " << spline.brent_bracket_itrs << std::endl;

  CHECK(mean_calls < 7.0);
  CHECK(mean_itrs < 4.0);
}
*/

