/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Revised Release on the IAPWS-95 formulation
  for the Thermodynamic Properties of Water and Steam.
*/

#include <array>

#include "Brent.hh"
#include "PowellHybrid.hh"

#include "IAPWS95.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* F(rho) = F3(rho) - p = 0
****************************************************************** */
struct Frho95 {
  Frho95(double p, double T, IAPWS95* eos) : p_(p), T_(T), eos_(eos), RT_(eos_->R * T_ / 1000.0) {};
  double operator()(double rho) const {
    double delta = rho / eos_->RHOC;

    double gd = eos_->ResidualPart(rho, T_)[1];
    double po = RT_ * (1.0 + delta * gd) * rho;
    return po - p_;
  }

  double p_, T_;
  IAPWS95* eos_;
  const double RT_;
};


/* ******************************************************************
* Calculate all properties for (p,T) input data
****************************************************************** */
std::tuple<Properties, Properties, Properties>
IAPWS95::ThermodynamicsPT(double p, double T)
{
  // initial guess and estimates of root brackets
  itrs_ = 20;
  double tol = 1e-9;
  double rho0 = eos97_.ThermodynamicsPT(p, T).rho; 
  double rhomin = rho0 * 0.995;
  double rhomax = rho0 * 1.005;

  Frho95 f(p, T, this);
  double rho = Utils::findRootBrent(f, rhomin, rhomax, tol, &itrs_);
  brent_root_itrs += itrs_;

  // refine soltution strategy by bracketing a root starting with a twice 
  // bigger bracket than before
  if (itrs_ < 0) {
    itrs_ = 20;
    auto [rhomin, rhomax] = Utils::bracketRootSymmetric(f, rho0, rho0 * 0.01, &itrs_);
    AMANZI_ASSERT(itrs_ >= 0);
    brent_bracket_itrs += itrs_;

    itrs_ = 20;
    rho = Utils::findRootBrent(f, rhomin, rhomax, tol, &itrs_);
    AMANZI_ASSERT(itrs_ > 0);
    brent_root_itrs += itrs_;
  }

  return ThermodynamicsRhoT(rho, T);
}


/* ******************************************************************
* Calculate all properties for (rho,T) input data
****************************************************************** */
std::tuple<Properties, Properties, Properties>
IAPWS95::ThermodynamicsRhoT(double rho, double T)
{
  Properties prop, liquid, vapor;
  prop = PopulateProperties(rho, T);

  // two-phase properties
  bool two_phase(false);
  double x(0.0), tol(1e-6);

  if (T < TC) {
    double rhol0, rhov0, rhol, rhov, p;
    rhol0 = DensityLiquid(T);
    rhov0 = DensityVapor(T);

    if (rhol0 > rho && rho > rhov0) {
      std::tie(rhol, rhov, p) = SaturationLineT(T, rhol0, rhov0);
      if (rhol * (1.0 + tol) > rho && rho > rhov * (1.0 - tol)) {
        liquid = PopulateProperties(rhol, T);
        vapor = PopulateProperties(rhov, T);

        liquid = ExtendProperties(rhol, liquid);
        vapor = ExtendProperties(rhov, vapor);

        // mean extensive properies
        double vl = 1.0 / rhol;
        double vv = 1.0 / rhov;
        x = (prop.v - vl) / (vv - vl);
        x = std::clamp(x, 0.0, 1.0);

        liquid.p = p;
        vapor.p = p;
        prop.p = p;

        prop.u = (1.0 - x) * liquid.u + x * vapor.u;
        prop.h = (1.0 - x) * liquid.h + x * vapor.h;
        prop.s = (1.0 - x) * liquid.s + x * vapor.s;
        prop.helmholtz = (1.0 - x) * liquid.helmholtz + x * vapor.helmholtz;

        two_phase = true;
      }
    } 

    if (!two_phase) {
      x = (rhol0 > rho) ? 0.0 : 1.0;
    }
  } else {
    x = 1.0;
  }

  prop = ExtendProperties(rho, prop);
  prop.x = x;

  if (!two_phase) {
    if (x == 0.0) liquid = prop;
    else vapor = prop;
  }

  return { prop, liquid, vapor };
}


/* ******************************************************************
* Calculate all properties for (p, h) input data
****************************************************************** */
std::tuple<Properties, Properties, Properties>
IAPWS95::ThermodynamicsPH(double p, double h)
{
  Properties prop, liquid, vapor;

  int phase = 1;
  if (p < PC) {
    const SaturationState& sat = SaturationLineP(p);

    if (h < sat.hl) {
      prop = PopulatePropertiesFromEntropy1(p, h);
      liquid = prop;
    } else if (h > sat.hv) {
      prop = PopulatePropertiesFromEntropy1(p, h);
      vapor = prop;
    } else {
      prop = PopulatePropertiesFromEntropy2(p, h, sat);
      liquid = prop;
      liquid.rho = sat.rhol;
      liquid.v = sat.vl;

      vapor = prop;
      vapor.rho = sat.rhov;
      vapor.v = sat.vv;
      phase = 2;
    }
  } else {
    prop = PopulatePropertiesFromEntropy1(p, h);
  }

  // extend state with kinematic properties 
  // for a homogeneous model, a default is logarithmic interpolation in volume fraction
  double T = prop.T;
  if (phase == 1) {
    double rho = prop.rho;
    prop.mu = Viscosity(rho, T);
    prop.k = ThermalConductivity(rho, T, prop);
  } else {
    double rhol = liquid.rho;
    double rhov = vapor.rho;
    double x = prop.x;
    double alpha = x / rhol / ((1 - x) / rhol + x / rhov);

    liquid.mu = Viscosity(rhol, T);
    vapor.mu = Viscosity(rhov, T);
    prop.mu = std::pow(liquid.mu, 1.0 - alpha) * std::pow(vapor.mu, alpha);

    liquid.k = ThermalConductivity(rhol, T, liquid);
    vapor.k = ThermalConductivity(rhov, T, vapor);
    prop.k = std::pow(liquid.k, 1.0 - alpha) * std::pow(vapor.k, alpha);

    prop.kt = std::pow(liquid.kt, 1.0 - alpha) * std::pow(vapor.kt, alpha);
  }

  return { prop, liquid, vapor };
}


/* ******************************************************************
* Populate state from the Helmholtz free energy
****************************************************************** */
Properties
IAPWS95::PopulateProperties(double rho, double T)
{
  Properties prop;

  const std::array<double, 6>& g0 = IdealGasPart(rho, T);
  const std::array<double, 6>& g = ResidualPart(rho, T);

  double delta = rho / RHOC;
  double tau = TC / T;

  double delta2 = delta * delta;
  double tau2 = tau * tau;
  double RT = R * T;

  double dg1 = delta * g[1];
  double dg3 = delta2 * g[3];
  double delta_tau_g4 = delta * tau * g[4];

  const double Z = 1.0 + dg1;
  const double D = 1.0 + 2.0 * dg1 + dg3;
  const double A = Z - delta_tau_g4;

  double g02 = g0[2] + g[2];
  double g05 = g0[5] + g[5];
  double tau_g02 = tau * g02;
  double tau2_g05 = tau2 * g05;

  prop.rho = rho;
  prop.T = T;

  prop.p = Z * RT * rho / 1000.0;
  prop.h = RT * (1.0 + tau_g02 + dg1);
  prop.u = RT * tau_g02;
  prop.s = R * (tau_g02 - g0[0] - g[0]);

  prop.cv = -R * tau2_g05;
  prop.cp = prop.cv + A * A / D;

  prop.ap = (1.0 - delta_tau_g4 / Z) / T;
  prop.av = (Z - delta_tau_g4) / (T * D);
  prop.bp = rho * (1.0 + (dg1 + dg3) / Z);

  prop.w = std::sqrt(1000.0 * RT * (D - A * A / tau2_g05));

  prop.v = 1.0 / rho;

  prop.helmholtz = RT * (g0[0] + g[0]);
  prop.gibbs = prop.helmholtz + 1000.0 * prop.p * prop.v;

  prop.kt = 1000.0 / (rho * RT * D);

  return prop;
}


/* ******************************************************************
* Populate state from entropy at its derivatives
****************************************************************** */
Properties
IAPWS95::PopulatePropertiesFromEntropy1(double p, double h)
{
  Properties prop;

  std::array<double, 6> a = EntropyDerivativesPH(p, h); 
  double s   = a[0];
  double sp  = a[1];
  double sh  = a[2];
  double spp = a[3];
  double sph = a[4];
  double shh = a[5];

  double T = 1.0 / sh;
  double Tp = -sph / (sh * sh);
  double Th = -shh / (sh * sh);

  double v = -sp / sh;
  double vp = -(spp * sh - sp * sph) / (sh * sh);
  double vh = -(sph * sh - sp * shh) / (sh * sh);

  double rhop = -vp / (v * v);
  double rhoh = -vh / (v * v);

  double drho_dp_s = rhop + v * rhoh;
  prop.cp = Th > 0.0 ? 1.0 / Th : std::numeric_limits<double>::infinity();
  prop.w = drho_dp_s > 0.0 ? std::sqrt(1.0 / drho_dp_s) : std::numeric_limits<double>::quiet_NaN();

  prop.p = p;
  prop.T = T;
  prop.u = h - p * v * 1000;
  prop.h = h;

  prop.rho = 1.0 / v;
}


/* ******************************************************************
* Populate state from entropy at its derivatives
****************************************************************** */
Properties
IAPWS95::PopulatePropertiesFromEntropy2(double p, double h, const SaturationState& sat)
{
  Properties prop;

  double dh = sat.hv - sat.hl;
  double dv = sat.vv - sat.vl;

  double x = (h - sat.hl) / dh;
  x = std::clamp(x, 0.0, 1.0);

  double dh_dp = sat.hv_p - sat.hl_p;
  double dv_dp = sat.vv_p - sat.vl_p;

  double xh = 1.0 / dh;
  double xp = -(sat.hl_p + x * dh_dp) / dh;

  double T = sat.Tsat;
  double Tp = sat.Tsat_p;
  double Th = 0.0;

  // volume and its derivatives
  double v = sat.vl + x * dv;
  double vp = sat.vl_p + x * dv_dp + xp * dv;
  double vh = xh * dv;

  double rhop = -vp / (v * v);
  double rhoh = -vh / (v * v);

  double s = sat.sl + x * (sat.sv - sat.sl);
  double sh = 1.0 / T;
  double sp = -v / T;

  double shh = 0.0;
  double sph = -Tp / (T * T);
  double spp = -vp / T + v * Tp / (T * T);

  prop.p = p;
  prop.T = T;
  prop.h = h;

  double drho_dp_s = rhop + v * rhoh;
  prop.w = drho_dp_s > 0.0 ? std::sqrt(1.0 / drho_dp_s) : 0.0;

  prop.rho = 1.0 / v;
  prop.v = v;
  prop.u = h - p * v * 1000;
  prop.s = s;
  prop.x = x;

  prop.helmholtz = prop.u - T * s;
  prop.gibbs = h - T * s;

  return prop;
}


/* ******************************************************************
* Ideal gas part of Helmholtz free energy
* http://www.iapws.org/relguide/IAPWS-95.html, Table 4
****************************************************************** */
std::array<double, 6>
IAPWS95::IdealGasPart(double rho, double T)
{
  double delta, tau;
  delta = rho / RHOC;
  tau = TC / T;

  double g, gd, gt, gdd, gdt(0.0), gtt;
  g = std::log(delta) + n0[1] + n0[2] * tau + n0[3] * std::log(tau);

  gd = 1.0 / delta;
  gdd = -gd / delta;
 
  gt = n0[2] + n0[3] / tau;
  gtt = -n0[3] / tau / tau;
  for (int i = 4; i < 9; ++i) {
    double tmp = 1.0 - std::exp(-gamma0[i - 4] * tau);
    g += n0[i] * std::log(tmp);
    gt += n0[i] * gamma0[i - 4] * (1.0 / tmp - 1.0);
    gtt -= n0[i] * gamma0[i - 4] * gamma0[i - 4] * (1.0 / tmp - 1.0) / tmp;
  } 

  return { g, gd, gt, gdd, gdt, gtt };
}


/* ******************************************************************
* Residual part of Hemholtz free energy
* http://www.iapws.org/relguide/IAPWS-95.html
****************************************************************** */
std::array<double, 6>
IAPWS95::ResidualPart(double rho, double T)
{
  double delta, tau;
  delta = rho / RHOC;
  tau = TC / T;

  double dpow[Nd_max];
  double tpow[Nt_max];
  double epow[Nc_max];

  dpow[0] = 1.0;
  for (int i = 1; i < Nd_max; ++i)
    dpow[i] = dpow[i - 1] * delta;

  tpow[0] = 1.0;
  for (int i = 1; i < Nt_max; ++i)
    tpow[i] = tpow[i - 1] * tau;

  epow[0] = 1.0;
  for (int i = 1; i < Nc_max; ++i)
    epow[i] = std::exp(-dpow[i]);

  double tpow1[7];
  tpow1[0] = std::pow(tau, -0.5);
  tpow1[1] = std::pow(tau, 0.875);
  tpow1[2] = tau;
  tpow1[3] = 1.0 / tpow1[0];
  tpow1[4] = std::pow(tau, 0.75);
  tpow1[5] = std::pow(tau, 0.375);
  tpow1[6] = tau;

  double g(0.0), gd(0.0), gt(0.0), gdd(0.0), gdt(0.0), gtt(0.0);

  // polynomial terms
  for (int i = 0; i < 7; ++i) {
    g += n1[i] * dpow[d1[i]] * tpow1[i];
    if (d1[i] > 0) gd += n1[i] * d1[i] * dpow[d1[i] - 1] * tpow1[i];
    if (d1[i] > 1) gdd += n1[i] * d1[i] * (d1[i] - 1) * dpow[d1[i] - 2] * tpow1[i];

    double tmp1 = dpow[d1[i]];
    double tmp2 = std::pow(tau, t1[i] - 1.0);
    gt += n1[i] * t1[i] * tmp1 * tmp2;
    if (t1[i] != 1.0) gtt += n1[i] * t1[i] * (t1[i] - 1.0) * tmp1 * std::pow(tau, t1[i] - 2.0);

    if (d1[i] > 0) gdt += n1[i] * d1[i] * t1[i] * dpow[d1[i] - 1] * tmp2;
  }

  // exponential terms
  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  for (int i = 0; i < 44; ++i) {
    int c = c2[i];
    int d = d2[i];
    int t = t2[i];
    double tmp1 = dpow[d] * epow[c];
    double tmp2 = tpow[t];

    g += n2[i] * tmp1 * tmp2;
    gd += n2[i] * tmp2 * dpow[d - 1] * (d - c * dpow[c]) * epow[c]; 
    if (d == 1) {
      gdd += n2[i] * tmp2 * epow[c] * c * dpow[c - 1] * (-1 + c * dpow[c] - c);
    } else {
      gdd += n2[i] * tmp2 * epow[c] * dpow[d - 2] * ((d - c * dpow[c]) * (d - 1 - c * dpow[c]) - c * c * dpow[c]);
    }

    gt += n2[i] * t * tmp1 * tpow[t - 1];
    if (t > 1) gtt += n2[i] * t * (t - 1) * tmp1 * tpow[t - 2];

    gdt += n2[i] * t * tpow[t - 1] * dpow[d - 1] * (d - c * dpow[c]) * epow[c]; 
  }

  // Gaussian terms
  double al, be, ga;
  for (int i = 0; i < 3; ++i) {
    int d = d3[i];
    int t = t3[i];
    al = alpha3[i];
    be = beta3[i];
    ga = gamma3[i];

    tmp1 = al * (delta - 1) * (delta - 1);
    tmp2 = be * (tau - ga) * (tau - ga);
    tmp3 = std::exp(-tmp1 - tmp2);
    tmp4 = dpow[d] * tpow[t];
    tmp5 = t / tau - 2 * be * (tau - ga);

    g += n3[i] * tmp4 * tmp3;

    // we know that d > 1, so the result is "symmetric" to gt and gtt
    gd += n3[i] * tmp3 * dpow[d - 1] * tpow[t] * (d - 2 * al * delta * (delta - 1));
    gdd += n3[i] * dpow[d - 2] * tpow[t] * tmp3 *
      (2 * al * delta * delta * (2 * tmp1 - 1.0) - 4 * d * al * delta * (delta - 1) + d * (d - 1));

    gt += n3[i] * tmp4 * tmp3 * (t / tau - 2 * be * (tau - ga));
    gtt += n3[i] * tmp4 * tmp3 * (tmp5 * tmp5 - t / tau / tau - 2 * be);

    gdt += n3[i] * tmp4 * tmp3 * (d / delta - 2 * al * (delta - 1)) * (t / tau - 2 * be * (tau - ga)); 
  }

  // other terms
  double theta;
  double del, deld(0.0), delt(0.0), deldd(0.0), deldt(0.0), deltt(0.0);
  double psi, psid, psit, psidd, psidt, psitt;

  for (int i = 0; i < 2; ++i) {
    tmp1 = (delta - 1) * (delta - 1);
    tmp2 = (tau - 1) * (tau - 1);

    tmp3 = std::pow(tmp1, 0.5 / beta4[i]);
    tmp4 = std::pow(tmp1, a4[i]);
    theta = (1.0 - tau) + A[i] * tmp3;
    del = theta * theta + B[i] * tmp4;
    psi = std::exp(-C[i] * tmp1 - D[i] * tmp2);

    psid = -2 * C[i] * (delta - 1) * psi;
    psit = -2 * D[i] * (tau - 1) * psi;

    psidd = 2 * C[i] * (2 * C[i] * tmp1 - 1) * psi;
    psitt = 2 * D[i] * (2 * D[i] * tmp2 - 1) * psi;
    psidt = 4 * C[i] * D[i] * (delta - 1) * (tau - 1) * psi;

    if (delta != 1.0) {
      tmp5 = theta * A[i] / beta4[i] * tmp3 / tmp1 + B[i] * a4[i] * tmp4 / tmp1;
      tmp6 = A[i] / beta4[i] * tmp3 / tmp1;
      deld = 2 * (delta - 1) * tmp5;
      deldd = 2 * tmp5 + tmp1 * (4 * B[i] * a4[i] * (a4[i] - 1) * tmp4 / tmp1 / tmp1 +
                                 2 * tmp6 * tmp6 + 
                                 4 * A[i] * theta / beta4[i] * (0.5 / beta4[i] - 1) * tmp3 / tmp1 / tmp1);

      deldt = -2 * tmp6 * (delta - 1);
    }

    delt = -2 * theta;
    deltt = 2.0;

    tmp3 = std::pow(del, b4[i]);
    tmp4 = b4[i] * std::pow(del, b4[i] - 1);
    tmp5 = b4[i] * (b4[i] - 1) * std::pow(del, b4[i] - 2);
    g += n4[i] * tmp3 * delta * psi;

    gd += n4[i] * (tmp3 * (psi + delta * psid) + tmp4 * deld * delta * psi);
    gdd += n4[i] * (tmp3 * (2 * psid + delta * psidd) + 
                    2 * tmp4 * deld * (psi + delta * psid) +
                    (tmp4 * deldd + tmp5 * deld * deld) * delta * psi);

    gt += n4[i] * delta * (psit * tmp3 + tmp4 * psi * delt);
    gtt += n4[i] * delta * (tmp3 * psitt + 2 * tmp4 * psit * delt + tmp5 * delt * delt * psi + tmp4 * deltt * psi);
    
    gdt += n4[i] * (tmp3 * (psit + delta * psidt) + delta * tmp4 * deld * psit + 
                    tmp4 * delt * (psi + delta * psid) + 
                    (tmp5 * deld * delt + tmp4 * deldt) * delta * psi);
  }

  residual_calls++;
  return { g, gd, gt, gdd, gdt, gtt };
}


/* ******************************************************************
* Derivatives are computed outside the saturation dome and inside 
* the metastable extension which requires direct call for computing 
* the ideal gas and residual parts of energy.
****************************************************************** */
FrhoT::Vector
FrhoT::operator()(FrhoT::Vector& x)
{
  const auto& g0 = eos_->IAPWS95::IdealGasPart(x[0], x[1]);
  const auto& gr = eos_->IAPWS95::ResidualPart(x[0], x[1]);

  double delta = x[0] / eos_->RHOC;
  double tau = eos_->TC / x[1];

  Vector r(x.size());
  r[0] = x[0] * eos_->R * x[1] * (1 + delta * gr[1]) / 1000 - p_;
  r[1] = eos_->R * x[1] * (1 + tau * (g0[2] + gr[2]) + delta * gr[1]) - h_;
  return r;
}


std::array<double, 6>
IAPWS95::EntropyDerivativesPH(double p, double h)
{
  auto [prop, liquid, vapor] = eos97_.ThermodynamicsPH(p, h); 
  double rho0 = prop.rho;
  double T0 = prop.T;

  // refine initial quess
  itrs_ = 10;
  double tol(1e-11);
  FrhoT f(p, h, this);
  FrhoT::Vector x0(2);
  x0[0] = rho0;
  x0[1] = T0;
  FrhoT::Vector sol = PowellHybrid(x0, f, &itrs_, tol);
  AMANZI_ASSERT(itrs_ >= 0);
  powell_root_itrs += itrs_;

  double rho = sol[0];
  double T = sol[1];
  return EntropyDerivativesPHbase(rho, T);
}


std::array<double, 6>
IAPWS95::EntropyDerivativesPHbase(double rho, double T)
{
  auto [prop, liquid, vapor] = ThermodynamicsRhoT(rho, T); 
  double s = prop.s;
  double cp = prop.cp;
  double alpha = prop.av;

  double sp = -1000.0 / (rho * T);  // p in MPa, h in kJ/kg
  double sh = 1.0 / T;

  double T2 = T * T;
  double rho2 = rho * rho;
  double Tp = 1000.0 * (alpha * T - 1.0) / (rho * cp);
  double rhop_h = rho2 / (prop.p * prop.bp) - rho * alpha * Tp;

  double shh = -1.0 / (T2 * cp);
  double sph = -Tp / T2;
  double spp = 1000.0 * (rhop_h / (rho2 * T) + Tp / (rho * T2));
  return {s, sp, sh, spp, sph, shh};
}


/* ******************************************************************
* Calculation of the two phase liquid and vapor boundaries.
****************************************************************** */
struct Frho2 {
  typedef Utils::VectorSTL Vector;

  Frho2(double T, IAPWS95* eos) : T_(T), eos_(eos) {};
  Vector operator()(Vector& x) {
    const auto& gl = eos_->IAPWS95::ResidualPart(x[0], T_);
    const auto& gv = eos_->IAPWS95::ResidualPart(x[1], T_);

    double delta_l = x[0] / eos_->RHOC;
    double delta_v = x[1] / eos_->RHOC;
    double jl = delta_l * (1.0 + delta_l * gl[1]);
    double jv = delta_v * (1.0 + delta_v * gv[1]);
    double kl = delta_l * gl[1] + gl[0] + std::log(delta_l);
    double kv = delta_v * gv[1] + gv[0] + std::log(delta_v);
  
    Vector r(x.size());
    r[0] = kv - kl;
    r[1] = jv - jl;
    return r;
  }
  double T_;
  IAPWS95* eos_;
};


std::tuple<double, double, double>
IAPWS95::SaturationLineT(double T, double rhol0, double rhov0)
{
  double psat, Tmin, tol(1e-11);
  Tmin = std::min(T, TC);

  itrs_ = 10;
  Frho2 f(T, this);
  Frho2::Vector x0(2);
  x0[0] = rhol0;
  x0[1] = rhov0;
  Frho2::Vector sol = PowellHybrid(x0, f, &itrs_, tol);
  powell_root_itrs += itrs_;

  if (sol[0] == sol[1]) {
    psat = PC;
  } else {
    double gl = IAPWS95::ResidualPart(sol[0], Tmin)[0];
    double gv = IAPWS95::ResidualPart(sol[1], Tmin)[0];

    psat = R * T * sol[0] * sol[1] / (sol[0] - sol[1]) * (gl - gv + std::log(sol[0] / sol[1])) / 1000.0;
  }
  return { sol[0], sol[1], psat };
}


/* ******************************************************************
* Calculation of the two phase liquid and vapor boundaries which
* returns more data.
****************************************************************** */
struct Frho3 {
  typedef Utils::VectorSTL Vector;

  Frho3(double p, IAPWS95* eos) : p_(p), eos_(eos) {};
  Vector operator()(Vector& x) {
    const auto& gl = eos_->IAPWS95::ResidualPart(x[0], x[2]);
    const auto& gv = eos_->IAPWS95::ResidualPart(x[1], x[2]);

    double delta_l = x[0] / eos_->RHOC;
    double delta_v = x[1] / eos_->RHOC;

    Vector r(x.size());
    r[0] = x[0] * eos_->R * x[2] * (1 + delta_l * gl[1]) / 1000 - p_;
    r[1] = x[1] * eos_->R * x[2] * (1 + delta_v * gv[1]) / 1000 - p_;
    r[2] = std::log(delta_l) + gl[0] + delta_l * gl[1] - (std::log(delta_v) + gv[0] + delta_v * gv[1]);
    return r;
  }
  double p_;
  IAPWS95* eos_;
};


SaturationState
IAPWS95::SaturationLineP(double p)
{
  SaturationState sat;

  if (p < PC) {
    double T0, T, rhol0, rhov0;
    T0 = eos97_.SaturationLineP(p);
    rhol0 = DensityLiquid(T0);
    rhov0 = DensityVapor(T0);

    itrs_ = 15;
    double tol(1e-11);
    Frho3 f(p, this);
    Frho3::Vector x0(3);
    x0[0] = rhol0;
    x0[1] = rhov0;
    x0[2] = T0;
    Frho3::Vector sol = PowellHybrid(x0, f, &itrs_, tol);
    AMANZI_ASSERT(itrs_ >= 0);
    powell_root_itrs += itrs_;

    sat.rhol = sol[0];
    sat.rhov = sol[1];
    T = sol[2];

    auto [prop1, liquid1, vapor1] = ThermodynamicsRhoT(sat.rhol, T);
    sat.hl = prop1.h;
    sat.sl = prop1.s;

    auto [prop2, liquid2, vapor2] = ThermodynamicsRhoT(sat.rhov, T);
    sat.hv = prop2.h;
    sat.sv = prop2.s;

    sat.vl = 1.0 / sat.rhol;
    sat.vv = 1.0 / sat.rhov;

    sat.p = p;
    sat.Tsat = T;
  }

  return sat;
}


/* ******************************************************************
* Auxiliary equation for the vapour pressure
* http://www.iapws.org/relguide/Supp-sat.html, Eq.1
****************************************************************** */
double
IAPWS95::VaporPressure(double T)
{
  if (T < TT) T = TT;
  else if (T > TC) T = TC;

  double tau = 1.0 - T / TC;
  double sum = 0.0;
  for (int i = 0; i < 6; ++i) {
    sum += PV_n[i] * std::pow(tau, PV_k[i]);
  }
  return std::exp(TC / T * sum) * PC;
}


/* ******************************************************************
* Equation for saturated liquid density, [kg/m3]
* http://www.iapws.org/relguide/Supp-sat.html, Eq.2
* http://www.iapws.org/relguide/Supp-sat.html, Eq.3
****************************************************************** */
double
IAPWS95::DensityLiquid(double T)
{
  if (T < TT) T = TT;
  else if (T > TC) T = TC;

  double tau = 1.0 - T / TC;
  double sum = 1.0;
  for (int i = 0; i < 6; ++i) {
    sum += RHOl_n[i] * std::pow(tau, RHOl_k[i] / 3);
  }
  return sum * RHOC;
}


double
IAPWS95::DensityVapor(double T)
{
  if (T < TT) T = TT;
  else if (T > TC) T = TC;

  double tau = 1.0 - T / TC;
  double sum = 0.0;
  for (int i = 0; i < 6; ++i) {
    sum += RHOv_n[i] * std::pow(tau, RHOv_k[i] / 6);
  }
  return std::exp(sum) * RHOC;
}


/* ******************************************************************
* Finalize properties
****************************************************************** */
Properties
IAPWS95::ExtendProperties(double rho, const Properties& prop_in)
{
  Properties prop;
  prop = prop_in;

  double T = prop.T;
  prop.mu = Viscosity(rho, T);
  prop.k = ThermalConductivity(rho, T, prop);

  return prop;
}


/* ******************************************************************
* i/o
****************************************************************** */
void
IAPWS95::Print(Properties& prop)
{
  std::cout << std::setprecision(12)
    << "============================" 
    << "\np = " << prop.p
    << "\nT = " << prop.T
    << "\nrho = " << prop.rho 
    << "\nv = " << prop.v
    << "\nh = " << prop.h << " = " << prop.u + prop.p * prop.v * 1000
    << "\nu = " << prop.u
    << "\ns = " << prop.s
    << "\nc_p = " << prop.cp
    << "\nc_v = " << prop.cv
    << "\nspeed of sound = " << prop.w
    << "\nalpha_v = " << prop.av
    << "\nalpha_p = " << prop.ap
    << "\nbeta_p = " << prop.bp
    << "\n\nHelmholtz = " << prop.helmholtz
    << "\nGibbs = " << prop.gibbs
    << "\n\nmu = " << prop.mu
    << "\nk = " << prop.k
    << "\nsigma = " << prop.sigma
    << "\nx = " << prop.x 
    << "\n=============================\n\n";
}
} // namespace AmanziEOS
} // namespace Amanzi
