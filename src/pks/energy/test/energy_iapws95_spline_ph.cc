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
  AmanziEOS::IAPWS95_RaggedSplinePH spline(plist);

  spline.InitializeSharedData();

  int count(0);
  double eps(2e-6);
  double dp(0.2), dh(5.0), hl, hr, der;

  spline.residual_calls = 0;
  const auto& opt = spline.GetOptions();
  for (double p = opt.p_min + dp; p < opt.p_max - dp; p += dp) {
    for (double h = opt.h_min + dh; h < opt.h_max - dh; h += dh) {
      const AmanziEOS::SaturationState& sat = eos95.SaturationLineP(p);
      if (spline.IsPhysical(p, h, sat)) {
        hl = h * (1.0 - eps);
        hr = h * (1.0 + eps);

        auto [prop1, liquid1, vapor1] = spline.ThermodynamicsPH(p, h);
        auto [prop3, liquid3, vapor3] = spline.ThermodynamicsPH(p, hr);
        auto [prop4, liquid4, vapor4] = spline.ThermodynamicsPH(p, hl);

        // Second-order entropy derivative with respect to h 
        der = (prop4.s - 2 * prop1.s + prop3.s) / (eps * eps * h * h);
        const auto& exact = spline.EntropyDerivativesPH(p, h);
        CHECK_CLOSE(der, exact[5], 2e-3 * std::max(1.0, std::fabs(exact[5])));

        count++;
      }
    }
  }

  auto calls = spline.residual_calls;
  auto itrs = spline.powell_root_itrs;
  double mean_calls = (double)calls / count;
  double mean_itrs = (double)itrs / count;
  std::cout << "Residual evaluations: " << calls << " or " << mean_calls << " times per sample\n";
  std::cout << "Powell's root find itrs: " << itrs << " or " << mean_itrs << " iterations per sample\n";
  std::cout << "Brent's bracket find itrs: " << spline.brent_bracket_itrs << std::endl;

  CHECK(mean_calls < 5.0);
  CHECK(mean_itrs < 5.0);
}

