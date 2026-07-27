#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>

#include "UnitTest++.h"
#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "IAPWS95.hh"
#include "IAPWS95_RaggedSpline.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziEOS;

void WriteHelmholtzPlotData(IAPWS95_RaggedSpline& spline,
                            IAPWS95& eos95,
                            double rho_min,
                            double rho_max,
                            double T_min,
                            double T_max,
                            int nrho,
                            int nT,
                            const std::string& filename)
{
  AMANZI_ASSERT(nrho >= 2 && nT >= 2);

  std::ofstream out(filename);
  out << std::setprecision(17);

  // Header
  out << "# rho T spline_phi exact_phi error\n";

  for (int j = 0; j < nT; ++j) {
    double T = T_min + (T_max - T_min) * (double)j / (nT - 1);

    for (int i = 0; i < nrho; ++i) {
      double rho = rho_min + (rho_max - rho_min) * (double)i / (nrho - 1);

      try {
        const auto& spline_value = spline.ResidualPart(rho, T);
        const auto& exact_value = eos95.ResidualPart(rho, T);
        double error = spline_value[0] - exact_value[0];

        out << rho << " " << T << " " << spline_value[0] << " " << exact_value[0] << " " << error << "\n";
      }
      catch (const std::exception&) {
        // Point is outside the ragged physical domain.
        // NaN preserves the rectangular plotting structure.
        const double nan = std::numeric_limits<double>::quiet_NaN();

        out << rho << " " << T << " " << nan << " " << nan << " " << nan << "\n";
      }
    }

    out << "\n";
  }
  std::cout << "Wrote plotting data to " << filename << '\n';
}


TEST(EOS_IAPWS95_SPLINE)
{
  IAPWS95_RaggedSpline::Options opt;
  opt.T_min = 300.0;
  opt.T_max = 900.0;
  opt.rho_min = 10.0;
  opt.rho_max = 1100.0;
  opt.initial_rho_intervals = 48;
  opt.initial_T_intervals = 48;
  opt.extension_cells = 4.0;
  opt.extension_weight = 0.1;

  Teuchos::ParameterList plist;
  IAPWS95_RaggedSpline spline(plist, opt);
  IAPWS95 eos95(plist);

  spline.CreateRaggedMesh();
  spline.BuildSplineCoefficients();

  // double rho(400.0), T(700.0);
  double rho(596.976), T(631.221);
  const auto& exact = eos95.ResidualPart(rho, T);
  const auto& approx = spline.ResidualPart(rho, T);

  for (int k = 0; k < 6; ++k) {
    std::cout << k << " exact = " << exact[k] << " spline = " << approx[k] << "\n";
    CHECK_CLOSE(exact[k], approx[k], 0.05);
  }

  // auto [prop, liquid, vapor] = eos95.ThermodynamicsRhoT(596.976, 631.221);
  // eos95.Print(prop); exit(0);

  const auto& mesh = spline.GetMesh();

  WriteHelmholtzPlotData(spline,
                         eos95,
                         mesh.rho_lines.front(),
                         mesh.rho_lines.back(),
                         mesh.T_lines.front(),
                         mesh.T_lines.back(),
                         301,
                         301,
                         "eos_iapws95_spline.dat");
}

