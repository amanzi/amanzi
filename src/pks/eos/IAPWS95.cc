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
  Frho95(double p, double T, IAPWS95* eos) : p_(p), T_(T), eos_(eos), RT_(IAPWS95::R * T / 1000.0) {};
  double operator()(double rho) const {
    double delta = rho / IAPWS95::RHOC;

    double gd = eos_->ResidualPartFirst(rho, T_)[1];
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
  int itrs = 20;
  double tol = 1e-9;
  double rho0 = std::get<0>(eos97_.ThermodynamicsPT(p, T)).rho; 
  double rhomin = rho0 * 0.995;
  double rhomax = rho0 * 1.005;

  Frho95 f(p, T, this);
  double rho = Utils::findRootBrent(f, rhomin, rhomax, tol, &itrs);
  brent_root_itrs += itrs;

  // refine soltution strategy by bracketing a root starting with a twice 
  // bigger bracket than before
  if (itrs < 0) {
    itrs = 20;
    auto [rhomin, rhomax] = Utils::bracketRootSymmetric(f, rho0, rho0 * 0.01, &itrs);
    AMANZI_ASSERT(itrs >= 0);
    brent_bracket_itrs += itrs;

    itrs = 20;
    rho = Utils::findRootBrent(f, rhomin, rhomax, tol, &itrs);
    AMANZI_ASSERT(itrs > 0);
    brent_root_itrs += itrs;
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
      prop.rgn = (int)Phase_t::CompressibleLiquid;
      liquid = prop;
    } else if (h > sat.hv) {
      prop = PopulatePropertiesFromEntropy1(p, h);
      prop.rgn = (int)Phase_t::Gas;
      vapor = prop;
    } else {
      prop = PopulatePropertiesFromEntropy2(p, h, sat);
      prop.rgn = (int)Phase_t::TwoPhases;
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
  prop.cp = prop.cv + R * A * A / D;

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

  // dh = T ds + 1000 v dp, because p is in MPa and h is in kJ/kg.
  double v = -sp / (1000.0 * sh);
  double vp = -(spp * sh - sp * sph) / (1000.0 * sh * sh);
  double vh = -(sph * sh - sp * shh) / (1000.0 * sh * sh);

  double rho = 1.0 / v;
  double rhop = -vp / (v * v);
  double rhoh = -vh / (v * v);

  double K = spp - sph * sph / shh;
  double rhop_T = 1000.0 * sh * K / (sp * sp);

  // Along isentrope ds = sp dp + sh dh = 0
  // hence (dh/dp)_s = -sp / sh = 1000 v
  double drho_dp_s = rhop + 1000.0 * v * rhoh;
  prop.cp = Th > 0.0 ? 1.0 / Th : std::numeric_limits<double>::infinity();
  prop.cv = -(sh * sh * spp - 2.0 * sp * sh * sph + sp * sp * shh) / (K * shh);
  prop.w = drho_dp_s > 0.0 ? std::sqrt(1.0e6 / drho_dp_s) : std::numeric_limits<double>::quiet_NaN();

  prop.p = p;
  prop.T = T;
  prop.u = h - p * v * 1000.0;
  prop.h = h;
  prop.s = s;

  prop.av = vh / (v * Th);
  prop.bp = rho * rho / (p * rhop_T);
  prop.kt = rhop_T / rho;

  prop.v = v;
  prop.rho = rho;
  return prop;
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
  double sp = -1000.0 * v / T;

  double shh = 0.0;
  double sph = -Tp / (T * T);
  double spp = 1000.0 * (-vp / T + v * Tp / (T * T));

  prop.p = p;
  prop.T = T;
  prop.h = h;

  double drho_dp_s = rhop + 1000.0 * v * rhoh;
  prop.w = drho_dp_s > 0.0 ? std::sqrt(1.0e6 / drho_dp_s) : 0.0;

  prop.rho = 1.0 / v;
  prop.v = v;
  prop.u = h - p * v * 1000.0;
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
* Derivatives are computed outside the saturation dome and inside 
* the metastable extension which requires direct call for computing 
* the ideal gas and residual parts of energy.
****************************************************************** */
FrhoT::Vector
FrhoT::operator()(FrhoT::Vector& x)
{
  const auto& g0 = eos_->IAPWS95::IdealGasPart(x[0], x[1]);
  const auto& gr = eos_->IAPWS95::ResidualPartFirst(x[0], x[1]);

  double delta = x[0] / IAPWS95::RHOC;
  double tau = IAPWS95::TC / x[1];

  Vector r(x.size());
  r[0] = x[0] * IAPWS95::R * x[1] * (1 + delta * gr[1]) / 1000 - p_;
  r[1] = IAPWS95::R * x[1] * (1 + tau * (g0[2] + gr[2]) + delta * gr[1]) - h_;
  return r;
}


std::array<double, 6>
IAPWS95::EntropyDerivativesPH(double p, double h)
{
  auto [prop, liquid, vapor] = eos97_.ThermodynamicsPH(p, h); 
  double rho0, T0;
  if (prop.rgn == 4) {
    double dhl = std::fabs(h - liquid.h);
    double dhv = std::fabs(h - vapor.h);
    if (dhl < dhv) {
      rho0 = liquid.rho;
      T0 = liquid.T;
    } else {
      rho0 = vapor.rho;
      T0 = vapor.T;
    }
  } else {
    rho0 = prop.rho;
    T0 = prop.T;
  }

  // refine initial quess
  int itrs = 20;
  double tol(1e-11);
  FrhoT f(p, h, this);
  FrhoT::Vector x0(2);
  x0[0] = rho0;
  x0[1] = T0;
  FrhoT::Vector sol = PowellHybrid(x0, f, &itrs, tol);
  AMANZI_ASSERT(itrs >= 0);
  powell_root_itrs += itrs;

  double rho = sol[0];
  double T = sol[1];
  return EntropyDerivativesPHbase(rho, T);
}


std::array<double, 6>
IAPWS95::EntropyDerivativesPHbase(double rho, double T)
{
  // use the homogeneous IAPWS95 state
  auto prop = PopulateProperties(rho, T);
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
    const auto& gl = eos_->IAPWS95::ResidualPartFirst(x[0], T_);
    const auto& gv = eos_->IAPWS95::ResidualPartFirst(x[1], T_);

    double delta_l = x[0] / IAPWS95::RHOC;
    double delta_v = x[1] / IAPWS95::RHOC;
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
  int itrs = 10;
  double psat, Tmin, tol(1e-11);
  Tmin = std::min(T, TC);

  Frho2 f(T, this);
  Frho2::Vector x0(2);
  x0[0] = rhol0;
  x0[1] = rhov0;
  Frho2::Vector sol = PowellHybrid(x0, f, &itrs, tol);
  AMANZI_ASSERT(itrs >= 0);
  powell_root_itrs += itrs;

  if (sol[0] == sol[1]) {
    psat = PC;
  } else {
    double gl = IAPWS95::ResidualPartFirst(sol[0], Tmin)[0];
    double gv = IAPWS95::ResidualPartFirst(sol[1], Tmin)[0];

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
    const auto& gl = eos_->IAPWS95::ResidualPartFirst(x[0], x[2]);
    const auto& gv = eos_->IAPWS95::ResidualPartFirst(x[1], x[2]);

    double delta_l = x[0] / IAPWS95::RHOC;
    double delta_v = x[1] / IAPWS95::RHOC;

    Vector r(x.size());
    r[0] = x[0] * IAPWS95::R * x[2] * (1 + delta_l * gl[1]) / 1000 - p_;
    r[1] = x[1] * IAPWS95::R * x[2] * (1 + delta_v * gv[1]) / 1000 - p_;
    r[2] = std::log(delta_l) + gl[0] + delta_l * gl[1] - (std::log(delta_v) + gv[0] + delta_v * gv[1]);
    return r;
  }
  double p_;
  IAPWS95* eos_;
};


SaturationState
IAPWS95::SaturationLineP(double p)
{
  SaturationState sat{};

  sat.p = p;

  if (p >= PC) {
    sat.Tsat = TC;
    sat.rhol = RHOC;
    sat.rhov = RHOC;
    sat.vl = 1.0 / RHOC;
    sat.vv = 1.0 / RHOC;
    sat.hl = HC;
    sat.hv = HC;

    auto prop = PopulateProperties(RHOC, TC);
    sat.sl = prop.s;
    sat.sv = prop.s;

    return sat;
  }

  if (p < PC) {
    double T0, T, rhol0, rhov0;
    T0 = eos97_.SaturationLineP(p);
    rhol0 = DensityLiquid(T0);
    rhov0 = DensityVapor(T0);

    int itrs = 15;
    double tol(1e-11);
    Frho3 f(p, this);
    Frho3::Vector x0(3);
    x0[0] = rhol0;
    x0[1] = rhov0;
    x0[2] = T0;
    Frho3::Vector sol = PowellHybrid(x0, f, &itrs, tol);
    AMANZI_ASSERT(itrs >= 0);
    powell_root_itrs += itrs;

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
