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
  Frho95(double p, double T, IAPWS95* eos) : p_(p), T_(T), eos_(eos) {};
  double operator()(double rho) const {
    double delta = rho / eos_->RHOC;

    double gd = eos_->ResidualPart(rho, T_)[1];
    double po = (1.0 + delta * gd) * eos_->R * T_ * rho / 1000.0;
    return po - p_;
  }

  double p_, T_;
  IAPWS95* eos_;
};


/* ******************************************************************
* Calculate all properties
****************************************************************** */
std::tuple<Properties, Properties, Properties>
IAPWS95::ThermodynamicsPT(double p, double T)
{
  itrs_ = 20;
  double tol = 1e-9;
  double rho0 = eos97_.ThermodynamicsPT(p, T).rho;
  double rhomin = rho0 * 0.905;
  double rhomax = rho0 * 1.005;

  Frho95 f(p, T, this);
  double rho = Utils::findRootBrent(f, rhomin, rhomax, tol, &itrs_);

  // refine soltution strategy starting with bracketing a root
  if (itrs_ < 0) {
    itrs_ = 20;
    auto [rhomin, rhomax] = Utils::bracketRoot(f, rho0, rho0 * 0.01, &itrs_);
    itrs_ = 20;
    rho = Utils::findRootBrent(f, rhomin, rhomax, tol, &itrs_);
  }

  return ThermodynamicsRhoT(rho, T);
}


std::tuple<Properties, Properties, Properties>
IAPWS95::ThermodynamicsRhoT(double rho, double T)
{
  const std::array<double, 6>& g0 = IdealGasPart(rho, T);
  const std::array<double, 6>& g = ResidualPart(rho, T);

  double delta = rho / RHOC;
  double tau = TC / T;

  double A = 1.0 + delta * g[1] - delta * tau * g[4];

  Properties prop, liquid, vapor;
  prop.rho = rho;
  prop.T = T;
  prop.p = (1 + delta * g[1]) * R * T * rho / 1000.0;
  prop.h = R * T * (1 + tau * (g0[2] + g[2]) + delta * g[1]);
  prop.u = R * T * tau * (g0[2] + g[2]);
  prop.s = R * (tau * (g0[2] + g[2]) - g0[0] - g[0]);
  prop.cv = -R * tau * tau * (g0[5] + g[5]);
  prop.cp = -R * tau * tau * (g0[5] + g[5]) + A * A / (1.0 + 2 * delta * g[1] + delta * delta * g[3]);
  prop.ap = (1 - delta * tau * g[4] / (1 + delta * g[1])) / T;
  prop.bp = rho * (1 + (delta * g[1] + delta * delta * g[3]) / (1 + delta * g[1]));
  prop.w = std::sqrt(R * T * 1000.0 * (1.0 + 2 * delta * g[1] + delta * delta * g[3]
                                           - A * A / (tau * tau * (g0[5] + g[5]))));

  prop.v = 1 / rho;
  prop.helmholtz = g[0];

  double x;
  if (T < TC) {
    double rhol, rhov, ps, p;
    rhol = DensityLiquid(T);
    rhov = DensityVapor(T);

    if (rhol > rho && rho > rhov) {
      std::tie(rhol, rhov, ps) = SaturationLine(T);
      if (rhol > rho && rho > rhov) {
        liquid = ExtendProperies(rhol, prop);
        vapor = ExtendProperies(rhov, prop);

        x = (1.0 / rho - 1.0 / rhol) / (1.0 / rhov - 1.0 / rhol);
        p = ps / 1000;

        liquid.p = p;
        vapor.p = p;
      }
    }
  } else {
    x = 1.0;
  }

  prop = ExtendProperies(rho, prop);
  prop.x = x;

  return { prop, liquid, vapor };
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

  return { g, gd, gt, gdd, gdt, gtt };
}


/* ******************************************************************
* Saturation calculation for two phase search FIXME
****************************************************************** */
struct Frho2 {
  typedef Utils::VectorSTL Vector;

  Frho2(double T, IAPWS95* eos) : T_(T), eos_(eos) {};
  Vector operator()(Vector& x) {
    const auto& gl = eos_->ResidualPart(x[0], T_);
    const auto& gv = eos_->ResidualPart(x[1], T_);

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
IAPWS95::SaturationLine(double T)
{
  double ps, Tmin;
  Tmin = std::min(T, TC);

  itrs_ = 10;
  Frho2 f(T, this);
  Frho2::Vector x0(2);
  x0[0] = DensityLiquid(T);
  x0[1] = DensityVapor(T);
  Frho2::Vector sol = PowellHybrid(x0, f, &itrs_);

  if (sol[0] == sol[1]) {
    ps = PC;
  } else {
    double gl = ResidualPart(sol[0], Tmin)[0];
    double gv = ResidualPart(sol[1], Tmin)[0];

    ps = R * T * sol[0] * sol[1] / (sol[0] - sol[1]) * (gl - gv + std::log(sol[0] / sol[1])) / 1000.0;
  }
  return { sol[0], sol[1], ps };
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
IAPWS95::ExtendProperies(double rho, const Properties& prop_in)
{
  Properties prop;
  prop = prop_in;

  double p = prop.p;
  double T = prop.T;

  prop.rho = rho;
  prop.v = 1.0 / rho;
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
    << "\nx = " << prop.x << "\n\n";
}
} // namespace AmanziEOS
} // namespace Amanzi
