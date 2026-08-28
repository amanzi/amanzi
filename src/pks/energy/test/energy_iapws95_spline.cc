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
#include "IAPWS95_RaggedSplineRhoT.hh"
#include "MeshFactory.hh"

TEST(HELMHOLTZ_IAPWS95_SPLINE_ENERGY_CONVEXITY)
{
  using namespace Amanzi;
  std::cout << "Test: Helmholtz energy spline properties" << std::endl;

  Teuchos::ParameterList plist;
  AmanziEOS::IAPWS95_RaggedSplineRhoT eos_spline(plist);

  eos_spline.InitializeSharedData();

  int n = 100;
  double eps(2e-6);
  double drho(2.2), dT(2.2), Tl, Tr, der;

  const auto& opt = eos_spline.GetOptions();
  for (double rho = opt.rho_min + drho; rho < opt.rho_max - drho; rho += drho) {
    for (double T = opt.T_min + dT; T < opt.T_max - dT; T += dT) {
      if (eos_spline.IsPhysical(rho, T)) {
        Tl = T * (1.0 - eps);
        Tr = T * (1.0 + eps);

        auto [prop1, liquid1, vapor1] = eos_spline.ThermodynamicsRhoT(rho, T);
        auto [prop3, liquid3, vapor3] = eos_spline.ThermodynamicsRhoT(rho, Tr);
        auto [prop4, liquid4, vapor4] = eos_spline.ThermodynamicsRhoT(rho, Tl);

        // Helmholtz energy is concave in T
        der = (prop4.helmholtz - 2 * prop1.helmholtz + prop3.helmholtz) / (eps * eps * T * T);
        CHECK(der < 0.0);
      }
    }
  }
}


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

