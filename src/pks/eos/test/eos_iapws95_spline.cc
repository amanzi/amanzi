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
#include "IAPWS95_RaggedSplineRhoT.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziEOS;

void WriteHelmholtzPlotData(IAPWS95_RaggedSplineRhoT& spline,
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

  std::array<double, 6> err_abs{}, err_rel{};
  for (int j = 0; j < nT; ++j) {
    double T = T_min + (T_max - T_min) * (double)j / (nT - 1);

    for (int i = 0; i < nrho; ++i) {
      double rho = rho_min + (rho_max - rho_min) * (double)i / (nrho - 1);

      if (spline.IsPhysical(rho, T)) {
        const auto& spline_value = spline.ResidualPart(rho, T);
        const auto& exact_value = eos95.ResidualPart(rho, T);
        double error = spline_value[0] - exact_value[0];

        for (int k = 0; k < 6; ++k) {
          double error = spline_value[k] - exact_value[k];
          err_abs[k] = std::max(err_abs[k], std::fabs(error));
          err_rel[k] = std::max(err_rel[k], std::fabs(error / std::max(1e-14, std::fabs(exact_value[k]))));
        }

        out << rho << " " << T << " " << spline_value[0] << " " << exact_value[0] << " " << error << "\n";
      } else {
        // Point is outside the ragged physical domain.
        // NaN preserves the rectangular plotting structure.
        const double nan = std::numeric_limits<double>::quiet_NaN();

        out << rho << " " << T << " " << nan << " " << nan << " " << nan << "\n";
      }
    }

    out << "\n";
  }
  std::cout << "\nWrote plotting data to " << filename << '\n';
  std::cout << "Spline resolutons: " << spline.GetMesh().x_lines.size() << " " << spline.GetMesh().y_lines.size() << std::endl;
  std::cout << "Error:  absolute     relative\n";
  for (int k = 0; k < 6; ++k) printf("%3d %12.8f %12.8f\n", k, err_abs[k], err_rel[k]); 
}


TEST(EOS_IAPWS95_SPLINE_RHO_T)
{
  IAPWS95_RaggedSplineRhoT::Options opt;
  opt.T_min = 280.0;
  opt.T_max = 950.0;
  opt.rho_min = 1.0;
  opt.rho_max = 1100.0;
  opt.initial_rho_intervals = 60;
  opt.initial_T_intervals = 70;
  opt.max_rho_intervals = 130;
  opt.max_T_intervals = 180;
  opt.extension_cells = 4.0;

  Teuchos::ParameterList plist;
  IAPWS95_RaggedSplineRhoT spline(plist, opt);
  IAPWS95 eos95(plist);

  spline.CreateRaggedMesh();
  const auto samples = spline.BuildSamples();
  spline.BuildSplineCoefficients(samples);

  // double rho(400.0), T(700.0);
  double rho(596.976), T(631.221);
  const auto& exact = eos95.ResidualPart(rho, T);
  const auto& approx = spline.ResidualPart(rho, T);

  for (int k = 0; k < 6; ++k) {
    printf("%3d exact=%12.8f spline=%12.8f\n", k, exact[k], approx[k]); 
    CHECK_CLOSE(exact[k], approx[k], 0.05);
  }

  // (dp/drho)_T = RT (1 + D) / 1000
  const auto& mesh = spline.GetMesh();
  double rho_min = mesh.x_lines.front();
  double rho_max = mesh.x_lines.back();
  double T_min = mesh.y_lines.front();
  double T_max = mesh.y_lines.back();

  std::ofstream out("field.dat");
  out << std::setprecision(17);

  int nT(301), nrho(301);
  const double nan = std::numeric_limits<double>::quiet_NaN();

  for (int j = 0; j < nT; ++j) {
    for (int i = 0; i < nrho; ++i) {
      double T = T_min + (T_max - T_min) * (double)j / (nT - 1);
      double rho = rho_min + (rho_max - rho_min) * (double)i / (nrho - 1);

      if (spline.IsPhysical(rho, T)) {
        double delta = rho / eos95.RHOC;
        double tau = eos95.TC / T;
        const auto& exact0 = eos95.IdealGasPart(rho, T);
        const auto& exactr = eos95.ResidualPart(rho, T);
        const auto& approx = spline.ResidualPart(rho, T);

        double Dex = 1.0 + 2 * delta * exactr[1] + delta * delta * exactr[3];
        double Dap = 1.0 + 2 * delta * approx[1] + delta * delta * approx[3];
        CHECK(Dap > 0.0);
        if (Dap < 0.0 || Dex < 0.0) std::cout << "rho=" << rho << "  T=" << T << "  D=" << Dex << " " << Dap << std::endl;
        // out << rho << " " << T << " " << Dap << " " << std::fabs(Dex - Dap) / std::max(Dex, Dap)<< "\n";

        double Nex = 1.0 + delta * exactr[1] - delta * tau * exactr[4];
        double Nap = 1.0 + delta * approx[1] - delta * tau * approx[4];
        // out << rho << " " << T << " " << Nap << " " << std::fabs(Nex - Nap) / std::max(Nex, Nap)<< "\n";

        double CPex = (-tau * tau * eos95.R * (exact0[5] + exactr[5]) + eos95.R * Nex * Nex / Dex) * 1000.0;
        double CPap = (-tau * tau * eos95.R * (exact0[5] + approx[5]) + eos95.R * Nap * Nap / Dap) * 1000.0;
        out << rho << " " << T << " " << CPap << " " << std::fabs(CPex - CPap) / std::max(CPex, CPap)<< "\n";
      }
      else {
        out << rho << " " << T << " " << nan << " " << nan << "\n";
      }
    }
  }
  out.close();

  WriteHelmholtzPlotData(spline,
                         eos95,
                         mesh.x_lines.front(),
                         mesh.x_lines.back(),
                         mesh.y_lines.front(),
                         mesh.y_lines.back(),
                         301,
                         301,
                         "eos_iapws95_spline.dat");
}

