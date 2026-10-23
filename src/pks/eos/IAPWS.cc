/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Implementation of transport properties for IAPWS ordinary water model.
*/

#include <algorithm>
#include <cmath>

#include "Brent.hh"
#include "PowellHybrid.hh"

#include "IAPWS97.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* Thermal conductivity
****************************************************************** */
double
IAPWS::ThermalConductivity(double rho, double T, Properties& prop)
{
  double k = ThermalConductivityBase(rho, T);
  if (k_enhancement_) k += ThermalConductivityCriticalEnhancement(prop);
  return k;
}


/* ******************************************************************
* http://www.iapws.org/relguide/ThCond.html, formulas (15)-(17)
****************************************************************** */
double
IAPWS::ThermalConductivityBase(double rho, double T) const
{
  double rhor, Tr;
  rhor = rho / RHOC;
  Tr = T / TC;

  // first factor
  const auto& n0 = thermal_cond_n0;

  double tmp(1.0), inva(1.0 / Tr);
  double k0(0.0), k1(0.0);
  for (int i = 0; i < N0_ThCond; ++i) {
    k0 += n0[i] * tmp;
    tmp *= inva;
  } 
  k0 = std::sqrt(Tr) / k0;
 
  // second factor
  double a(inva - 1.0), b(rhor - 1.0); 
  double apow[N0_ThCond];
  double bpow[N1_ThCond];

  apow[0] = 1.0;
  for (int i = 1; i < N0_ThCond; ++i)
    apow[i] = apow[i - 1] * a;

  bpow[0] = 1.0;
  for (int i = 1; i < N1_ThCond; ++i)
    bpow[i] = bpow[i - 1] * b;

  const auto& n1 = thermal_cond_n1;

  for (int i = 0; i < N0_ThCond; ++i) {
    for (int j = 0; j < N1_ThCond; ++j) {
      k1 += n1[i][j] * apow[i] * bpow[j];
    }
  }
  k1 = std::exp(rhor * k1);
  return 1e-3 * k0 * k1;
}


/* ******************************************************************
* Formula (18)
****************************************************************** */
double
IAPWS::ThermalConductivityCriticalEnhancement(const Properties& prop) const
{
  // This cutoff is NOT part of the IAPWS formulation.
  // if (std::fabs(prop.T - Tc) > 100.0 && std::fabs(prop.rho - rhoc) > 200.0) return 0.0;

  constexpr double Tbar_ref = 1.5;
  constexpr double R95 = 0.46151805;  // Release uses gas constant from IAPWS95 formulation
  constexpr double crit_max = 1.0e+13;  // saveguard maximum value for zeta and cp

  // Reduced temperature and density.
  double Tr = prop.T / TC;
  double rhor = prop.rho / RHOC;

  // Dimensionless isothermal density derivative: zeta = (d delta / d pi)_T
  constexpr double scale = PC / RHOC;
  double zeta = scale * prop.rho * prop.kt;
  double zeta_ref = ThermalConductivityZetaFunction_(rhor);
  if (zeta < 0.0 || zeta > crit_max) zeta = crit_max;

  // Critical susceptibility increment:
  // IAPWS explicitly requires negative Delta_chi to be replaced by zero.
  double delta_chi = rhor * (zeta - zeta_ref * Tbar_ref / Tr);
  if (delta_chi <= 0.0) return 0.0;

  // Correlation length
  constexpr double nu_over_gamma = nu / gamma;
  double xi = xi0_k * std::pow(delta_chi / Gamma0, nu_over_gamma);
  double y = xi / qD_inv_k;

  // IAPWS specifies Z(y) = 0 below this value to avoid numerical truncation problems
  if (y < 1.2e-7) return 0.0;

  // Crossover function Z(y)
  double kappa_inv = prop.cv / prop.cp;
  double A = 1.0 / y + (y * y) / (3.0 * rhor * rhor);
  double omega = (1.0 - kappa_inv) * std::atan(y) + kappa_inv * y;
  double omega0 = -std::expm1(-1.0 / A);  // for better accuracy when argument is small
  double Z = 2.0 / (M_PI * y) * (omega - omega0);

  // Small negative values are possible only from roundoff
  if (Z <= 0.0) return 0.0;

  double cp_bar = prop.cp / R95;
  double mu_bar = prop.mu / 1.0e-6;
  if (cp_bar < 0.0 || cp_bar > crit_max) cp_bar = crit_max;

  return 1e-3 * Lambda * rhor * cp_bar * Tr * Z / mu_bar;
}


/* ******************************************************************
* Formula (26)
****************************************************************** */
double
IAPWS::ThermalConductivityZetaFunction_(double rhor) const
{
  for (int i = 0; i < N2_ThCond; ++i) {
    if (rhor <= thermal_cond_rho[i]) {
      double zeta = thermal_cond_zeta[M2_ThCond - 1][i];
      for (int j = M2_ThCond - 2; j >= 0; --j) {
        zeta = zeta * rhor + thermal_cond_zeta[j][i];
      }
      return 1.0 / zeta;
    }
  }
}


/* ******************************************************************
* Surface tension
* http://www.iapws.org/relguide/Surf-H2O.html
* Validity range 248.15 <= T <= TC
****************************************************************** */
double
IAPWS::SurfaceTension(double T)
{
  if (248.15 <= T && T <= TC) {
    double tau = 1.0 - T / TC;
    return 235.8e-3 * std::pow(tau, 1.256) * (1.0 - 0.625 * tau);
  }
  return -1.0;
}


/* ******************************************************************
* Dynamic viscosity
****************************************************************** */
double
IAPWS::Viscosity(double rho, double T, const Properties& prop)
{
  double mu = ViscosityBase(rho, T);
  if (mu_enhancement_) mu *= ViscosityCriticalEnhancement(prop);
  return mu;
}


/* ******************************************************************
* http://www.iapws.org/relguide/viscosity.html, formulas (10)-(12)
****************************************************************** */
double
IAPWS::ViscosityBase(double rho, double T) const
{
  double rhor = rho / RHOC;
  double Tr = T / TC;

  // first factor
  const auto& n0 = viscosity_n0;

  double tmp(1.0), inva(1.0 / Tr);
  double mu0(0.0), mu1(0.0);
  for (int i = 0; i < N0_Visc; ++i) {
    mu0 += n0[i] * tmp;
    tmp *= inva;
  } 
  mu0 = 100.0 * std::sqrt(Tr) / mu0;
 
  // second factor
  const auto& k = viscosity_k1;
  const auto& l = viscosity_l1;
  const auto& n1 = viscosity_n1;

  double a(inva - 1.0), b(rhor - 1.0); 
  double apow[N1_Visc_kmax];
  double bpow[N1_Visc_lmax];

  apow[0] = 1.0;
  for (int i = 1; i < N1_Visc_kmax; ++i)
    apow[i] = apow[i - 1] * a;

  bpow[0] = 1.0;
  for (int i = 1; i < N1_Visc_lmax; ++i)
    bpow[i] = bpow[i - 1] * b;

  for (int i = 0; i < N1_Visc; ++i) {
    mu1 += n1[i] * apow[k[i]] * bpow[l[i]];
  }
  mu1 = std::exp(rhor * mu1);

  return 1e-6 * mu0 * mu1;
}


/* ******************************************************************
* Derivatives
****************************************************************** */
double
IAPWS::ViscosityBaseDerivativeRho(double rho, double T) const
{
  double rhor = rho / RHOC;
  double Tr = T / TC;

  // first factor depends on T only
  const auto& n0 = viscosity_n0;

  double tmp(1.0), inva(1.0 / Tr);
  double mu0(0.0), mu1, A(0.0), B(0.0);
  for (int i = 0; i < N0_Visc; ++i) {
    mu0 += n0[i] * tmp;
    tmp *= inva;
  } 
  mu0 = 100.0 * std::sqrt(Tr) / mu0;
 
  // second factor
  const auto& k = viscosity_k1;
  const auto& l = viscosity_l1;
  const auto& n1 = viscosity_n1;

  double a(inva - 1.0), b(rhor - 1.0); 
  double apow[N1_Visc_kmax];
  double bpow[N1_Visc_lmax];

  apow[0] = 1.0;
  for (int i = 1; i < N1_Visc_kmax; ++i)
    apow[i] = apow[i - 1] * a;

  bpow[0] = 1.0;
  for (int i = 1; i < N1_Visc_lmax; ++i)
    bpow[i] = bpow[i - 1] * b;

  for (int i = 0; i < N1_Visc; ++i) {
    double tmp = n1[i] * apow[k[i]];
    A += tmp * bpow[l[i]];
    if (l[i] > 0) B += tmp * l[i] * bpow[l[i] - 1];
  }
  mu1 = std::exp(rhor * A);
  double dmu1_drho = mu1 * (A + rhor * B) / RHOC;

  return 1e-6 * mu0 * dmu1_drho;
}


double
IAPWS::ViscosityBaseDerivativeT(double rho, double T) const
{
  double rhor = rho / RHOC;
  double Tr = T / TC;

  // first factor
  const auto& n0 = viscosity_n0;

  double tmp(1.0), inva(1.0 / Tr);
  double mu0, mu1, A(0.0), B(0.0);
  for (int i = 0; i < N0_Visc; ++i) {
    A += n0[i] * tmp;
    tmp *= inva;
    if (i > 0) B += n0[i] * i * tmp;
  } 
  double C = std::sqrt(Tr);
  mu0 = 100.0 * C / A;
  double dmu0_dT = 100.0 * (0.5 / (C * A) + C * B / (A * A));
 
  // second factor
  const auto& k = viscosity_k1;
  const auto& l = viscosity_l1;
  const auto& n1 = viscosity_n1;

  double a(inva - 1.0), b(rhor - 1.0); 
  double apow[N1_Visc_kmax];
  double bpow[N1_Visc_lmax];

  apow[0] = 1.0;
  for (int i = 1; i < N1_Visc_kmax; ++i)
    apow[i] = apow[i - 1] * a;

  bpow[0] = 1.0;
  for (int i = 1; i < N1_Visc_lmax; ++i)
    bpow[i] = bpow[i - 1] * b;

  A = 0.0;
  B = 0.0;
  for (int i = 0; i < N1_Visc; ++i) {
    double tmp = n1[i] * bpow[l[i]];
    A += tmp * apow[k[i]];
    if (k[i] > 0) B += tmp * k[i] * apow[k[i] - 1];
  }
  mu1 = std::exp(rhor * A);
  double dmu1_dT = -mu1 * rhor * B / (Tr * Tr);

  return 1e-6 * (mu0 * dmu1_dT + dmu0_dT * mu1) / TC;
}


/* ******************************************************************
* IAPWS R12-08 viscosity critical multiplier mu2_bar.
* Use R15-11 zeta_ref polynomial as an optional speed approximation
****************************************************************** */
double
IAPWS::ViscosityCriticalEnhancement(const Properties& prop) const
{
  // optional performace heuristic
  // if (prop.T < 640.0 || prop.T > 655.0 || prop.rho < 200.0 || prop.rho > 450.0) return 1.0;

  constexpr double Tbar_ref = 1.5;

  double Tr = prop.T / TC;
  double rhor = prop.rho / RHOC;

  // dimensionless isothermal susceptibility
  double zeta = (PC / RHOC) * prop.rho * prop.kt;
  double zeta_ref = ThermalConductivityZetaFunction_(rhor);
  double delta_chi = rhor * (zeta - zeta_ref * Tbar_ref / Tr);
  if (delta_chi <= 0.0) return 1.0;

  // formula (20)
  double xi = xi0_mu * std::pow(delta_chi / Gamma0, nu / gamma);
  double qCxi = xi / qC_inv_mu;
  double qDxi = xi / qD_inv_mu;

  double Y = 0.0;
  if (xi <= xi_switch) {
    // formula (15), small-correlation-length branch
    double qDxi2 = qDxi * qDxi;
    double qDxi5 = qDxi2 * qDxi2 * qDxi;
    Y = 0.2 * qCxi * qDxi5 * (1.0 - qCxi + qCxi * qCxi - (765.0 / 504.0) * qDxi2);
  } else {
    // formulas (16)-(19), large-correlation-length branch
    double qCxi2 = qCxi * qCxi;
    double qCxi3 = qCxi2 * qCxi;
    double psi_D = std::acos( 1.0 / std::sqrt(1.0 + qDxi * qDxi));

    double ratio = std::abs((qCxi - 1.0) / (qCxi + 1.0));
    double w = std::sqrt(ratio) * std::tan(0.5 * psi_D);

    double L;
    if (qCxi > 1.0) { // For this branch 0 <= w < 1.
      L = std::log((1.0 + w) / (1.0 - w));
    } else {
      L = 2.0 * std::atan(std::abs(w));
    }

    double abs_term = std::pow(std::fabs(qCxi2 - 1.0), 1.5);
    Y = (1.0 / 12.0) * std::sin(3.0 * psi_D)
      - (1.0 / (4.0 * qCxi)) * std::sin(2.0 * psi_D)
      + (1.0 / qCxi2) * (1.0 - 1.25 * qCxi2) * std::sin(psi_D)
      - (1.0 / qCxi3) * ((1.0 - 1.5 * qCxi2) * psi_D - abs_term * L);
  }

  return std::exp(x_mu * Y);
}

} // namespace AmanziEOS
} // namespace Amanzi

