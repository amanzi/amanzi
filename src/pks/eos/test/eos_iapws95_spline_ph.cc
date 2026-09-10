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
#include "IAPWS95_RaggedSplinePH.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziEOS;

double Dfun(const std::array<double, 6>& der, double R) {
  double sp = der[1];
  double spp = der[3];
  double sph = der[4];
  double shh = der[5];
  return sp * sp / (IAPWS95::R * (spp - sph * sph / shh));
}


/* ******************************************************************
* 
****************************************************************** */
void WriteEntropyPlotData(IAPWS95_RaggedSplinePH& spline,
                          IAPWS95& eos95,
                          double p_min,
                          double p_max,
                          double h_min,
                          double h_max,
                          int np,
                          int nh,
                          const std::string& filename)
{
  std::ofstream out(filename);
  out << std::setprecision(17);

  // Header
  out << "# p h spline_s exact_s error\n";

  int count(0);
  std::array<double, 6> err_l2{}, err_abs{}, err_rel{}, err{};
  std::array<double, 6> err_p_max{}, err_h_max{};

  for (int j = 0; j < np; ++j) {
    double p = p_min + (p_max - p_min) * (double)j / (np - 1);

    for (int i = 0; i < nh; ++i) {
      double h = h_min + (h_max - h_min) * (double)i / (nh - 1);

      const auto& sat = eos95.SaturationLineP(p);
      if (spline.IsPhysical(p, h, sat)) {
        const auto& spline_value = spline.EntropyDerivativesPH(p, h);
        const auto& exact_value = eos95.EntropyDerivativesPH(p, h);

        double Dex = Dfun(exact_value, IAPWS95::R);
        double Dap = Dfun(spline_value, IAPWS95::R);
        // CHECK(Dap > 0.0);
        if (Dap < 0.0 || Dex < 0.0) count++;
        out << p << " " << h << " " << Dap << " " << Dex << " " << Dap - Dex << "\n";

        for (int k = 0; k < 6; ++k) {
          double tmp = spline_value[k] - exact_value[k];
          if (std::fabs(tmp) > err_abs[k]) {
            err_p_max[k] = p;
            err_h_max[k] = h;
          }
          err[k] = tmp;
          err_l2[k] += err[k] * err[k];
          err_abs[k] = std::max(err_abs[k], std::fabs(err[k]));
          err_rel[k] = std::max(err_rel[k], std::fabs(err[k] / std::max(1e-14, std::fabs(exact_value[k]))));
        }
        // out << p << " " << h << " " << spline_value[3] << " " << exact_value[3] << " " << err[3] << "\n";
      }
      else {
        // Point is outside the ragged physical domain.
        // NaN preserves the rectangular plotting structure.
        const double nan = std::numeric_limits<double>::quiet_NaN();
        out << p << " " << h << " " << nan << " " << nan << " " << nan << "\n";
      }
    }

    out << "\n";
  }

  std::cout << "Wrote plotting data to " << filename << '\n';
  std::cout << "Spline resolutions: " << spline.GetMesh().x_lines.size() << " " << spline.GetMesh().y_lines.size() << std::endl;
  std::cout << "Number of negative D: " << count << std::endl;
  CHECK(count == 0);
  std::cout << "Error:  mean       max        relative_max   p_max      h_max\n";
  for (int k = 0; k < 6; ++k) {
    err_l2[k] = std::sqrt(err_l2[k] / np / nh);
    printf("%3d %12.8f %12.8f %12.8f %12.8f %12.8f\n", k, err_l2[k], err_abs[k], err_rel[k], err_p_max[k], err_h_max[k]); 
    if (k != 3) CHECK(err_l2[k] < 0.01);
  }
}


/* ******************************************************************
* Graphics support
****************************************************************** */
void WriteSaturationPlotData(IAPWS95_RaggedSplinePH& spline,
                             IAPWS95& eos95,
                             const std::string& filename)
{
  std::ofstream out(filename);
  out << std::setprecision(17);

  const auto& mesh = spline.GetMesh();
  const auto& lines = mesh.x_lines;
  for (int n = 0; n < lines.size(); ++n) {
    double p = lines[n];
    if (p <= eos95.PC) {
      const auto& sat = eos95.SaturationLineP(p);
      out << p << " " << sat.hl << std::endl;
    }
  }
  for (int n = lines.size() - 1; n >=0; --n) {
    double p = lines[n];
    if (p <= eos95.PC) {
      const auto& sat = eos95.SaturationLineP(p);
      out << p << " " << sat.hv << std::endl;
    }
  }
  out.close();
}


/* ******************************************************************
* Graphics support
****************************************************************** */
struct Frho95t {
  Frho95t(double p, double rho, IAPWS95* eos) : p_(p), rho_(rho), eos_(eos) {};
  double operator()(double T) const {
    double delta = rho_ / eos_->RHOC;

    double gd = eos_->ResidualPart(rho_, T)[1];
    double po = (1.0 + delta * gd) * eos_->R * T * rho_ / 1000.0;
    return po - p_;
  }

  double p_, rho_;
  IAPWS95* eos_;
};


double FindLiquidMetastableBoundary(IAPWS95_RaggedSplinePH& spline,
                                    IAPWS95* eos95,
                                    double p,
                                    const SaturationState& sat)
{
  double rho_hi = sat.rhol;
  double T_hi = sat.Tsat;
  double D_hi = spline.StabilityFactor(rho_hi, T_hi);

  if (!(D_hi > 0.0)) return sat.hl;

  // Follow the constant-pressure branch until D = (dp/drho)_T changes sign.
  double rho_lo = rho_hi;
  double T_lo = T_hi;
  double D_lo = D_hi;

  bool bracketed = false;

  for (int k = 0; k < 1000; ++k) {
    double rho1 = rho_hi * 0.995;

    // 1. Solve p(rho1, T1) = p staying on the continuation branch by centering
    //    the temperature bracket at the previous value T_hi.
    int itrs = 20;
    Frho95t f(p, rho1, eos95);

    auto [Tmin, Tmax] = Utils::bracketRootSymmetric(f, T_hi, T_hi * 0.01, &itrs);
    if (itrs < 0) break;

    itrs = 50;
    double T1 = Utils::findRootBrent(f, Tmin, Tmax, 1.0e-10, &itrs);
    if (itrs < 0) break;

    double D1 = spline.StabilityFactor(rho1, T1);
    if (!std::isfinite(D1)) break;

    // Stable endpoint: D_hi > 0
    // Unstable endpoint: D_lo <= 0
    if (D1 <= 0.0) {
      rho_lo = rho1;
      T_lo = T1;
      D_lo = D1;

      bracketed = true;
      break;
    }

    rho_hi = rho1;
    T_hi = T1;
    D_hi = D1;
  }

  if (!bracketed) return sat.hl;

  // 2. Refine the spinodal: p(rho,T) = p, D(rho,T) = 0 by bisecting in rho.
  //    For every rho, T is obtained from the constant-pressure equation.
  double rho_tol = 1.0e-10 * std::max(1.0, sat.rhol);
  double D_scale = std::max(std::abs(D_hi), 1.0);
  double D_tol = 1.0e-10 * D_scale;

  double rho_spin = 0.5 * (rho_lo + rho_hi);
  double T_spin = 0.5 * (T_lo + T_hi);

  for (int k = 0; k < 80; ++k) {
    rho_spin = 0.5 * (rho_lo + rho_hi);

    // Interpolate T between the two continuation states.
    // This gives an excellent center for the inner pressure solve.
    double theta = (rho_spin - rho_lo) / (rho_hi - rho_lo);
    double Tguess = T_lo + theta * (T_hi - T_lo);

    int itrs = 20;
    Frho95t f(p, rho_spin, eos95);

    auto [Tmin, Tmax] = Utils::bracketRootSymmetric(f, Tguess, std::max(0.1, 0.002 * Tguess), &itrs);
    if (itrs < 0) return sat.hl;

    itrs = 50;
    T_spin = Utils::findRootBrent(f, Tmin, Tmax, 1.0e-11, &itrs);
    if (itrs < 0) return sat.hl;

    double D_spin = spline.StabilityFactor(rho_spin, T_spin);
    if (!std::isfinite(D_spin)) return sat.hl;

    if (std::abs(D_spin) <= D_tol || rho_hi - rho_lo <= rho_tol) {
      auto prop = eos95->PopulateProperties(rho_spin, T_spin);
      return prop.h;
    }

    if (D_spin > 0.0) {
      // Still on stable/metastable side.
      rho_hi = rho_spin;
      T_hi = T_spin;
      D_hi = D_spin;
    } else {
      // Crossed the spinodal.
      rho_lo = rho_spin;
      T_lo = T_spin;
      D_lo = D_spin;
    }
  }

  auto prop = eos95->PopulateProperties(rho_spin, T_spin);
  return prop.h;
}


double FindVaporMetastableBoundary(IAPWS95_RaggedSplinePH& spline,
                                   IAPWS95* eos95,
                                   double p,
                                   const SaturationState& sat)
{
  double rho_lo = sat.rhov;
  double T_lo = sat.Tsat;
  double D_lo = spline.StabilityFactor(rho_lo, T_lo);

  if (!(D_lo > 0.0)) return sat.hv;

  double rho_hi = rho_lo;
  double T_hi = T_lo;
  double D_hi = D_lo;

  bool bracketed = false;

  // Find a density interval where D changes sign.
  for (int k = 0; k < 1000; ++k) {
    double rho1 = rho_lo * 1.005;

    int itrs = 20;
    Frho95t f(p, rho1, eos95);

    auto [Tmin, Tmax] = Utils::bracketRootSymmetric(f, T_lo, T_lo * 0.01, &itrs);
    if (itrs < 0) break;

    itrs = 50;
    double T1 = Utils::findRootBrent(f, Tmin, Tmax, 1.0e-10, &itrs);
    if (itrs < 0) break;

    double D1 = spline.StabilityFactor(rho1, T1);
    if (!std::isfinite(D1)) break;

    if (D1 <= 0.0) {
      rho_hi = rho1;
      T_hi = T1;
      D_hi = D1;
      bracketed = true;
      break;
    }

    rho_lo = rho1;
    T_lo = T1;
    D_lo = D1;
  }

  if (!bracketed) return sat.hv;

  double rho_tol = 1.0e-10 * std::max(1.0, sat.rhov);
  double D_scale = std::max(std::abs(D_lo), 1.0);
  double D_tol = 1.0e-10 * D_scale;

  double rho_spin = 0.5 * (rho_lo + rho_hi);
  double T_spin = 0.5 * (T_lo + T_hi);

  for (int k = 0; k < 80; ++k) {
    rho_spin = 0.5 * (rho_lo + rho_hi);
    double theta = (rho_spin - rho_lo) / (rho_hi - rho_lo);
    double Tguess = T_lo + theta * (T_hi - T_lo);

    int itrs = 20;
    Frho95t f(p, rho_spin, eos95);

    auto [Tmin, Tmax] = Utils::bracketRootSymmetric(f, Tguess, std::max(0.1, 0.002 * Tguess), &itrs);
    if (itrs < 0) return sat.hv;

    itrs = 50;
    T_spin = Utils::findRootBrent(f, Tmin, Tmax, 1.0e-11, &itrs);
    if (itrs < 0) return sat.hv;

    double D_spin = spline.StabilityFactor(rho_spin, T_spin);
    if (!std::isfinite(D_spin)) return sat.hv;

    if (std::abs(D_spin) <= D_tol || rho_hi - rho_lo <= rho_tol) {
      auto prop = eos95->PopulateProperties(rho_spin, T_spin);
      return prop.h;
    }

    if (D_spin > 0.0) {
      // Still on vapor-side metastable branch.
      rho_lo = rho_spin;
      T_lo = T_spin;
      D_lo = D_spin;
    } else {
      // Crossed vapor spinodal.
      rho_hi = rho_spin;
      T_hi = T_spin;
      D_hi = D_spin;
    }
  }

  auto prop = eos95->PopulateProperties(rho_spin, T_spin);
  return prop.h;
}


void WriteSpinodalPlotData(IAPWS95_RaggedSplinePH& spline,
                           IAPWS95& eos95,
                           const std::string& filename)
{
  std::ofstream out(filename);
  out << std::setprecision(17);

  const auto& mesh = spline.GetMesh();
  const auto& lines = mesh.x_lines;
  for (int n = 0; n < lines.size(); ++n) {
    double p = lines[n];
    if (p <= eos95.PC) {
      const auto& sat = eos95.SaturationLineP(p);
      double ext_hl = FindLiquidMetastableBoundary(spline, &eos95, p, sat);
      out << p << " " << ext_hl << std::endl;
    }
  }
  for (int n = lines.size() - 1; n >=0; --n) {
    double p = lines[n];
    if (p <= eos95.PC) {
      const auto& sat = eos95.SaturationLineP(p);
      double ext_hv = FindVaporMetastableBoundary(spline, &eos95, p, sat);
      out << p << " " << ext_hv << std::endl;
    }
  }
  out.close();
}


/* ******************************************************************
* 
****************************************************************** */
void WriteSamplesData(const std::vector<IAPWS95_RaggedSplinePH::Sample>& samples,
                      const std::string& filename)
{
  std::ofstream out(filename);
  out << std::setprecision(17);
  for (int n = 0; n < samples.size(); ++n) {
    out << samples[n].p << " " << samples[n].h << " " << samples[n].rho << std::endl;
  }
  out.close();
}


/* ******************************************************************
* 
****************************************************************** */
TEST(EOS_IAPWS95_SPLINE_P_H)
{
  IAPWS95_RaggedSplinePH::Options opt{};
  opt.p_min = 0.2;
  opt.p_max = 50.0;
  opt.h_min = 500.0;
  opt.h_max = 3600.0;
  opt.initial_p_intervals = 60;
  opt.initial_h_intervals = 70;
  opt.max_p_intervals = 150;
  opt.max_h_intervals = 180;
  opt.metastable_stability_fraction = 0.5;

  Teuchos::ParameterList plist;
  IAPWS95_RaggedSplinePH spline(plist, opt);
  IAPWS95 eos95(plist);

  const auto& samples = spline.InitializeSharedData();

  std::ofstream out("field.dat");
  out << std::setprecision(17);

  int np(301), nh(301);
  const double nan = std::numeric_limits<double>::quiet_NaN();

  for (int j = 0; j < np; ++j) {
    for (int i = 0; i < nh; ++i) {
      double p = opt.p_min + (opt.p_max - opt.p_min) * (double)j / (np - 1);
      double h = opt.h_min + (opt.h_max - opt.h_min) * (double)i / (nh - 1);

      if (spline.IsPhysical(p, h)) {
        const auto& [prop0, liquid0, vapor0] = eos95.ThermodynamicsPH(p, h);
        const auto& [prop1, liquid1, vapor1] = spline.ThermodynamicsPH(p, h);
        out << p << " " << h << " " << prop0.w << " " << prop1.w << "\n";
      }
      else {
        out << p << " " << h << " " << nan << " " << nan << "\n";
      }
    }
  }
  out.close();


  const auto& mesh = spline.GetMesh();
  WriteEntropyPlotData(spline,
                       eos95,
                       mesh.x_lines.front(),
                       mesh.x_lines.back(),
                       mesh.y_lines.front(),
                       mesh.y_lines.back(),
                       301,
                       301,
                       "eos_iapws95_spline_ph.dat");

  WriteSaturationPlotData(spline, eos95, "phase_boundary.dat");
  WriteSpinodalPlotData(spline, eos95, "spinodal_boundary.dat");
  WriteSamplesData(samples, "samples.dat");
}

