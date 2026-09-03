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

#include "Brent.hh"
#include "dbc.hh"
#include "IAPWS95.hh"
#include "IAPWS95_RaggedSplineRhoT.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziEOS;

/* ******************************************************************
* Graphics support
****************************************************************** */
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

  std::array<double, 6> err_l2{}, err_abs{}, err_rel{}, err{};
  for (int j = 0; j < nT; ++j) {
    double T = T_min + (T_max - T_min) * (double)j / (nT - 1);

    for (int i = 0; i < nrho; ++i) {
      double rho = rho_min + (rho_max - rho_min) * (double)i / (nrho - 1);

      if (spline.IsPhysical(rho, T)) {
        const auto& spline_value = spline.ResidualPart(rho, T);
        const auto& exact_value = eos95.ResidualPart(rho, T);
        std::array<double, 6> error;

        for (int k = 0; k < 6; ++k) {
          err[k] = spline_value[k] - exact_value[k];
          err_l2[k] += err[k] * err[k];
          err_abs[k] = std::max(err_abs[k], std::fabs(err[k]));
          err_rel[k] = std::max(err_rel[k], std::fabs(err[k] / std::max(1e-14, std::fabs(exact_value[k]))));
        }

        out << rho << " " << T << " " << spline_value[0] << " " << exact_value[0] << " " << err[0] << "\n";
      } else {
        // Point is outside the ragged physical domain.
        // NaN preserves the rectangular plotting structure.
        const double nan = std::numeric_limits<double>::quiet_NaN();

        out << rho << " " << T << " " << nan << " " << nan << " " << nan << "\n";
      }
    }

    out << "\n";
  }
  out.close();

  std::cout << "\nWrote plotting data to " << filename << '\n';
  std::cout << "Spline resolutions: " << spline.GetMesh().x_lines.size() << " " << spline.GetMesh().y_lines.size() << std::endl;
  std::cout << "Error:  mean       absolute     relative\n";
  for (int k = 0; k < 6; ++k) {
    err_l2[k] = std::sqrt(err_l2[k] / nrho / nT);
    printf("%3d %12.8f %12.8f %12.8f\n", k, err_l2[k], err_abs[k], err_rel[k]); 
  }
}


/* ******************************************************************
* Graphics support
****************************************************************** */
void WriteSaturationPlotData(IAPWS95_RaggedSplineRhoT& spline, const std::string& filename)
{
  std::ofstream out(filename);
  out << std::setprecision(17);

  const auto& saturation = spline.GetSaturation();
  for (int n = 0; n < saturation.size(); ++n) {
    out << saturation[n].rho_l << " " << saturation[n].T << std::endl;
  }
  for (int n = saturation.size() - 1; n >=0; --n) {
    out << saturation[n].rho_v << " " << saturation[n].T << std::endl;
  }
  out.close();
}


/* ******************************************************************
* Graphics support
****************************************************************** */
void WriteSpinodalPlotData(IAPWS95_RaggedSplineRhoT& spline, const std::string& filename)
{
  std::ofstream out(filename);
  out << std::setprecision(17);

  const auto& saturation = spline.GetSaturation();
  for (int n = 0; n < saturation.size(); ++n) {
    double T = saturation[n].T;
    double rho_hi = saturation[n].rho_l;
    double D_hi = spline.StabilityFactor(rho_hi, T);

    double rho_lo = rho_hi;
    double rho_min = saturation[n].rho_v;

    while (rho_lo > rho_min) {
      rho_lo *= 0.99;
      double D_lo = spline.StabilityFactor(rho_lo, T);

      if (D_lo <= 0.0) {
        int itrs = 20;
        double rho = Utils::findRootBrent([&](double rho) {
                                            return spline.StabilityFactor(rho, T);
                                          },
                                          rho_lo,
                                          rho_hi,
                                          1e-10,
                                          &itrs);
        AMANZI_ASSERT(itrs >= 0);
        out << rho << " " << saturation[n].T << std::endl;
        break;
      }

      rho_hi = rho_lo;
      D_hi = D_lo;
    }
  }

  for (int n = saturation.size() - 1; n >=0; --n) {
    double T = saturation[n].T;
    double rho_lo = saturation[n].rho_v;
    double D_lo = spline.StabilityFactor(rho_lo, T);

    double rho_hi = rho_lo;
    double rho_max = saturation[n].rho_l;

    while (rho_hi < rho_max) {
      rho_hi *= 1.01;
      double D_hi = spline.StabilityFactor(rho_hi, T);

      if (D_hi <= 0.0) {
        int itrs = 20;
        double rho = Utils::findRootBrent([&](double rho) {
                                            return spline.StabilityFactor(rho, T);
                                          },
                                          rho_lo,
                                          rho_hi,
                                          1e-10,
                                          &itrs);
        AMANZI_ASSERT(itrs >= 0);
        out << rho << " " << saturation[n].T << std::endl;
        break;
      }

      rho_lo = rho_hi;
      D_lo = D_hi;
    }
  }
  out.close();
}


/* ******************************************************************
* 
****************************************************************** */
void WriteSamplesData(const std::vector<IAPWS95_RaggedSplineRhoT::Sample>& samples,
                      const std::string& filename)
{
  std::ofstream out(filename);
  out << std::setprecision(17);
  for (int n = 0; n < samples.size(); ++n) {
    out << samples[n].rho << " " << samples[n].T << std::endl;
  }
  out.close();
}


/* ******************************************************************
* 
****************************************************************** */
TEST(EOS_IAPWS95_SPLINE_RHO_T)
{
  IAPWS95_RaggedSplineRhoT::Options opt;
  opt.T_min = 280.0;
  opt.T_max = 950.0;
  opt.rho_min = 0.5;
  opt.rho_max = 1100.0;
  opt.initial_rho_intervals = 60;
  opt.initial_T_intervals = 70;
  opt.max_rho_intervals = 130;
  opt.max_T_intervals = 180;
  opt.extension_cells = 4.0;
  opt.anisotropic_refinement_fraction = 0.01;

  Teuchos::ParameterList plist;
  {
    /*
    IAPWS95_RaggedSplineRhoT spline1(plist, opt);
    IAPWS95_RaggedSplineRhoT spline2(plist, opt);

    spline1.InitializeSharedData();
    spline2.InitializeSharedData();
    CHECK(std::addressof(spline1.GetMesh()) == std::addressof(spline2.GetMesh()));
    return; 
    */
  }

  IAPWS95_RaggedSplineRhoT spline(plist, opt);
  IAPWS95 eos95(plist);

  const auto samples = spline.InitializeSharedData();

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
        out << rho << " " << T << " " << Dap << " " << std::fabs(Dex - Dap) / std::max(Dex, Dap)<< "\n";

        double Nex = 1.0 + delta * exactr[1] - delta * tau * exactr[4];
        double Nap = 1.0 + delta * approx[1] - delta * tau * approx[4];
        // out << rho << " " << T << " " << Nap << " " << std::fabs(Nex - Nap) / std::max(Nex, Nap)<< "\n";

        double CPex = (-tau * tau * eos95.R * (exact0[5] + exactr[5]) + eos95.R * Nex * Nex / Dex) * 1000.0;
        double CPap = (-tau * tau * eos95.R * (exact0[5] + approx[5]) + eos95.R * Nap * Nap / Dap) * 1000.0;
        // out << rho << " " << T << " " << CPap << " " << std::fabs(CPex - CPap) / std::max(CPex, CPap)<< "\n";
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

  WriteSaturationPlotData(spline, "phase_boundary.dat");
  WriteSpinodalPlotData(spline, "spinodal_boundary.dat");
  WriteSamplesData(samples, "samples.dat");
}

