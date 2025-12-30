/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Revised Release on the IAPWS Industrial Formulation 1997
  for the Thermodynamic Properties of Water and Steam.
*/

#include "Brent.hh"
#include "PowellHybrid.hh"

#include "IAPWS97.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* F(rho) = F3(rho) - p = 0
****************************************************************** */
struct Frho {
  Frho(double p, double T, IAPWS97* eos) : p_(p), T_(T), eos_(eos) {};
  double operator()(double rho) const { return eos_->Region3(rho, T_).p - p_; }

  double p_, T_;
  IAPWS97* eos_;
};

struct Fh {
  Fh(double p, double h, int id, IAPWS97* eos) : p_(p), h_(h), id_(id), eos_(eos) {};
  double operator()(double T) const {
    if (id_ == 1) {
      return eos_->Region1(p_, T).h - h_;
    } else if (id_ == 2) {
      return eos_->Region2(p_, T).h - h_;
    } else if (id_ == 5) {
      return eos_->Region5(p_, T).h - h_;
    }
    return 0.0;
  }

  int id_;
  double p_, h_;
  IAPWS97* eos_;
};

struct FhT {
  typedef Utils::VectorSTL Vector;

  FhT(double p, double h, IAPWS97* eos) : p_(p), h_(h), eos_(eos) {};
  Vector operator()(Vector& x) {
    Vector r(x.size());
    r[0] = eos_->Region3(x[0], x[1]).h - h_;
    r[1] = eos_->Region3(x[0], x[1]).p - p_;
    return r;
  }

  double p_, h_;
  IAPWS97* eos_;
};


/* ******************************************************************
* Calculate all properties
****************************************************************** */
Properties
IAPWS97::ThermodynamicsPT(double p, double T)
{
  Properties prop;
  prop.p = p;
  prop.T = T;

  int rgn = RegionIdPT(p, T);
  itrs_ = 0;

  switch (rgn) {
    case 1: {
      prop = Region1(p, T);
      break;
    }
    case 2: {
      prop = Region2(p, T);
      break;
    }
    case 3: {
      double rho;
      if (T == TC && p == PC) {
        rho = RHOC;
      } else {
        double v = Region3_BackwardPT(p, T);
        double rho0 = 1.0 / v;

        itrs_ = 10;
        double tol = 1e-8;
        double rhomin = rho0 * 0.999;
        double rhomax = rho0 * 1.001;
        Frho f(p, T, this);
        rho = Utils::findRootBrent(f, rhomin, rhomax, tol, &itrs_);
      }
      prop = Region3(rho, T);
      break;
    }
    case 5: {
      prop = Region5(p, T);
      break;
    }
    default: {
      Exceptions::amanzi_throw("Invalid P/T range");
    }
  }

  prop.rho = 1.0 / prop.v;
  prop.phase = PhaseId(p, T, prop.rgn, prop.x);

  return prop;
}


/* ******************************************************************
* Calculate all properties
****************************************************************** */
std::tuple<Properties, Properties, Properties>
IAPWS97::ThermodynamicsPH(double p, double h)
{
  Properties prop;
  prop.p = p;
  prop.h = h;

  int rgn = RegionIdPH(p, h);
  prop.rgn = rgn;

  double tol(1e-8), T, T0, Tmin, Tmax;
  itrs_ = 0;

  switch (rgn) {
    case 1: {
      T0 = Region1_BackwardPH(p, h);

      itrs_ = 10;
      Tmin = T0 * 0.999;
      Tmax = T0 * 1.001;
      Fh f(p, h, 1, this);
      T = Utils::findRootBrent(f, Tmin, Tmax, tol, &itrs_);
      prop = Region1(p, T);
      break;
    }
    case 2: {
      T0 = Region2_BackwardPH(p, h);

      itrs_ = 10;
      Tmin = T0 * 0.999;
      Tmax = T0 * 1.001;
      Fh f(p, h, 2, this);
      T = Utils::findRootBrent(f, Tmin, Tmax, tol, &itrs_);
      prop = Region2(p, T);
      break;
    }
    case 3: {
      double v0 = Region3_BackwardPH_v(p, h);
      double T0 = Region3_BackwardPH_T(p, h);

      itrs_ = 10;
      FhT f(p, h, this);
      FhT::Vector x0(2);
      x0[0] = 1.0 / v0;
      x0[1] = T0;
      FhT::Vector sol = PowellHybrid(x0, f, &itrs_);
      prop = Region3(sol[0], sol[1]);
      break;
    }
    case 4: {
      double h1, h2, x;
      T = SaturationLineP(p);
      if (T <= 623.15) {
        h1 = Region1(p, T).h;
        h2 = Region2(p, T).h;
      } else {
        h1 = Region4(p, 0.0).h;
        h2 = Region4(p, 1.0).h;
      }
      x = (h - h1) / (h2 - h1);
      prop = Region4(p, x);
      break;
    }
    case 5: {
      itrs_ = 10;
      Tmin = 1500.0;
      Tmax = 2273.15;
      Fh f(p, h, 5, this);
      T = Utils::findRootBrent(f, Tmin, Tmax, tol, &itrs_);
      prop = Region5(p, T);
      break;
    }
    default: {
      Exceptions::amanzi_throw("Invalid P/H range");
    }
  }

  Properties liquid, gas;

  // only liquid phase
  if (prop.x == 0.0) { 
    liquid = ExtendProperies(prop);
    liquid.sigma = SurfaceTension(prop.T);
  }
  // only gas phase
  else if (prop.x == 1.0) {
    gas = ExtendProperies(prop);
  }
  // two phases
  else {
    double T = prop.T;
    Properties prop_l, prop_g;

    if (623.15 < T && T <= TC) {
      double rhol = 1.0 / Region3_BackwardPX(p, 0.0, T);  
      prop_l = Region3(rhol, T);

      double rhov = 1.0 / Region3_BackwardPX(p, 1.0, T);
      prop_g = Region3(rhov, T);
    } else {
      prop_l = Region1(p, T);
      prop_g = Region2(p, T);
    }
    prop_l.x = prop.x;
    prop_l.rgn = prop.rgn;
    prop_g.x = prop.x;
    prop_g.rgn = prop.rgn;

    liquid = ExtendProperies(prop_l);
    gas = ExtendProperies(prop_g);

    liquid.sigma = SurfaceTension(T);
  }

  return { prop, liquid, gas };
}


/* ******************************************************************
* Region 1, Section 5, formula 7
****************************************************************** */
Properties
IAPWS97::Region1(double p, double T)
{
  Properties prop;

  double pr = p / 16.53;
  double Tr = 1386 / T;

  double a = 7.1 - pr;
  double b = Tr - 1.222;
  double invb = 1.0 / b;

  double apow[N1_kmax];
  double bpow[N1_lmax];
  double bpow_inv[N1_lmax_inv];

  apow[0] = 1.0;
  for (int i = 1; i < N1_kmax; ++i)
    apow[i] = apow[i - 1] * a;

  bpow[0] = 1.0;
  for (int i = 1; i < N1_lmax; ++i)
    bpow[i] = bpow[i - 1] * b;

  bpow_inv[0] = 1.0;
  for (int i = 1; i < N1_lmax_inv; ++i)
    bpow_inv[i] = bpow_inv[i - 1] * invb;

  const auto& k = rgn1k;
  const auto& l = rgn1l;
  const auto& n = rgn1n;

  double g(0.0), gp(0.0), gpp(0.0), gt(0.0), gtt(0.0), gpt(0.0), tmp;
  for (int i = 0; i < N1; ++i) {
    tmp = (l[i] >= 0) ? bpow[l[i]] : bpow_inv[-l[i]];
    g += n[i] * apow[k[i]] * tmp;
    if (k[i] > 0) gp -= n[i] * k[i] * apow[k[i] - 1] * tmp;
    if (k[i] > 1) gpp += n[i] * k[i] * (k[i] - 1) * apow[k[i] - 2] * tmp;

    if (l[i] != 0) {
      tmp = (l[i] > 0) ? bpow[l[i] - 1] : bpow_inv[-l[i] + 1];
      gt += n[i] * l[i] * apow[k[i]] * tmp;
      if (k[i] > 0) gpt -= n[i] * k[i] * l[i] * apow[k[i] - 1] * tmp;
    }
    if (l[i] != 0 && l[i] != 1) {
      tmp = (l[i] > 0) ? bpow[l[i] - 2] : bpow_inv[-l[i] + 2];
      gtt += n[i] * l[i] * (l[i] - 1) * apow[k[i]] * tmp;
    }
  }

  double A = (gp - Tr * gpt) * (gp - Tr * gpt);

  prop.p = p;
  prop.T = T;
  prop.v = (pr / p) * gp * R * T / 1000;
  prop.h = (Tr * T) * gt * R;
  prop.u = R * T * (Tr * gt - pr * gp);
  prop.s = R * (Tr * gt - g);
  prop.cp = -R * Tr * Tr * gtt;
  prop.cv = R * (-Tr * Tr * gtt + A / gpp);
  prop.w = std::sqrt(R * T * 1000 * gp * gp / (A / (Tr * Tr * gtt) - gpp));
  prop.kt = -(pr / p) * gpp / gp;
  prop.av = (1 - Tr * gpt / gp) / T;
  prop.x = 0.0;
  prop.rgn = 1;

  // derived properties
  prop.rho = 1.0 / prop.v;

  return prop;
}


/* ******************************************************************
* Backward map for Region 1: (p, h) -> T, R7-97 (2012), formula 11.
****************************************************************** */
double
IAPWS97::Region1_BackwardPH(double p, double h)
{
  const auto& k = inv1k_ph;
  const auto& l = inv1l_ph;
  const auto& n = inv1n_ph;

  double pr = p;
  double hr = h / 2500.0;

  double a = pr;
  double b = hr + 1.0;

  double apow[N1ph_kmax];
  double bpow[N1ph_lmax];

  apow[0] = 1.0;
  for (int i = 1; i < N1ph_kmax; ++i) {
    apow[i] = apow[i - 1] * a;
  }

  bpow[0] = 1.0;
  for (int i = 1; i < N1ph_lmax; ++i) {
    bpow[i] = bpow[i - 1] * b;
  }

  double T(0.0);
  for (int i = 0; i < N1ph; ++i) {
    T += n[i] * apow[k[i]] * bpow[l[i]];
  }
  return T;
}


/* ******************************************************************
* Region 2, Section 6, formula 15
****************************************************************** */
Properties
IAPWS97::Region2(double p, double T)
{
  Properties prop;

  double pr = p;
  double Tr = 540.0 / T;
  double a = pr;
  double b = Tr;

  // ideal gas part of the dimensionless Gibbs free energy
  const auto& no = rgn2n_o;

  double go, gop, gopp, got, gott, gopt(0.0);
  go = std::log(pr);
  gop = 1.0 / pr;
  gopp = -gop / pr;
  
  double b2, b3, bneg1, bneg2, bneg3, bneg4, bneg5, bneg6, bneg7;
  b2 = b * b;
  b3 = b2 * b;
  bneg1 = 1.0 / b;
  bneg2 = bneg1 * bneg1;
  bneg3 = bneg2 * bneg1;
  bneg4 = bneg3 * bneg1;
  bneg5 = bneg4 * bneg1;
  bneg6 = bneg5 * bneg1;
  bneg7 = bneg6 * bneg1;

  go += no[0] + no[1] * b + no[2] * bneg5 + no[3] * bneg4 + no[4] * bneg3
      + no[5] * bneg2 + no[6] * bneg1 + no[7] * b2 + no[8] * b3;

  got = no[1] - 5 * no[2] * bneg6 - 4 * no[3] * bneg5 - 3 * no[4] * bneg4
      - 2 * no[5] * bneg3 - no[6] * bneg2 + 2 * no[7] * b + 3 * no[8] * b2;

  gott = 30 * no[2] * bneg7 + 20 * no[3] * bneg6 + 12 * no[4] * bneg5
      + 6 * no[5] * bneg4 + 2 * no[6] * bneg3 + 2 * no[7] + 6 * no[8] * b;

  // residual part of the Gibbs free energy
  a = pr;
  b = Tr - 0.5;

  double apow[N2r_kmax];
  double bpow[N2r_lmax];

  apow[0] = 1.0;
  for (int i = 1; i < N2r_kmax; ++i)
    apow[i] = apow[i - 1] * a;

  bpow[0] = 1.0;
  for (int i = 1; i < N2r_lmax; ++i)
    bpow[i] = bpow[i - 1] * b;
 
  const auto& k = rgn2k_r;
  const auto& l = rgn2l_r;
  const auto& nr = rgn2n_r;

  double gr(0.0), grp(0.0), grpp(0.0), grt(0.0), grtt(0.0), grpt(0.0), tmp;
  for (int i = 0; i < N2r; ++i) {
    tmp = bpow[l[i]];
    gr += nr[i] * apow[k[i]] * tmp;
    if (k[i] > 0) grp += nr[i] * k[i] * apow[k[i] - 1] * tmp;
    if (k[i] > 1) grpp += nr[i] * k[i] * (k[i] - 1) * apow[k[i] - 2] * tmp;

    tmp = apow[k[i]];
    if (l[i] > 0) grt += nr[i] * l[i] * tmp * bpow[l[i] - 1];
    if (l[i] > 1) grtt += nr[i] * l[i] * (l[i] - 1) * tmp * bpow[l[i] - 2];

    if (l[i] > 0 && k[i] > 0) grpt += nr[i] * k[i] * l[i] * apow[k[i] - 1] * bpow[l[i] - 1];
  }

  double A = (1 + pr * grp - Tr * pr * grpt) * (1 + pr * grp - Tr * pr * grpt);
  double B = 1 + 2 * pr * grp + pr * pr * grp * grp;

  prop.p = p;
  prop.T = T;
  prop.v = (pr / p) * (gop + grp) * R * T / 1000;
  prop.h = Tr * (got + grt) * R * T;
  prop.u = R * T * (Tr * (got + grt) - pr * (gop + grp));
  prop.s = R * (Tr * (got + grt) - (go + gr));
  prop.cp = -R * Tr * Tr * (gott + grtt);
  prop.cv = R * (-Tr * Tr * (gott + grtt) - A / (1 - pr * pr * grpp));
  prop.w = std::sqrt(R * T * 1000 * B / (1 - pr * pr * grpp + A / Tr / Tr / (gott + grtt)));
  prop.kt = (1 - pr * pr * grpp) / (1 + pr * grp) / p;
  prop.av = (1 + pr * grp - Tr * pr * grpt) / (1 + pr * grp) / T;
  prop.x = 1.0;
  prop.rgn = 2;

  // derived properties
  prop.rho = 1.0 / prop.v;

  return prop;
}


/* ******************************************************************
* Backward map for Region 2: (p, h) -> T, R7-97 (2012), formula 22-24
****************************************************************** */
double
IAPWS97::Region2_BackwardPH(double p, double h)
{
  double T, Tsat, hf;

  if (p <= 4.0) {
    T = Region2_SubregionFormulaA(p, h);
  } else if (4.0 < p && p <= 6.546699678) {
    T = Region2_SubregionFormulaB(p, h);
  } else {
    hf = BoundaryLine2b2c(p);
    if (h >= hf) {
      T = Region2_SubregionFormulaB(p, h);
    } else {
      T = Region2_SubregionFormulaC(p, h);
    }
  }

  if (p <= 22.064) {
    Tsat = SaturationLineP(p);
    T = std::max(T, Tsat);
  }
  return T;
}


/* ******************************************************************
* Backward map for subregion 2a: formula 22.
****************************************************************** */
double
IAPWS97::Region2_SubregionFormulaA(double p, double h)
{
  const auto& k = inv2k_a;
  const auto& l = inv2l_a;
  const auto& n = inv2n_a;

  double pr = p;
  double hr = h / 2000.0;

  double a = pr;
  double b = hr - 2.1;

  double apow[N2a_kmax];
  double bpow[N2a_lmax];

  apow[0] = 1.0;
  for (int i = 1; i < N2a_kmax; ++i) {
    apow[i] = apow[i - 1] * a;
  }

  bpow[0] = 1.0;
  for (int i = 1; i < N2a_lmax; ++i) {
    bpow[i] = bpow[i - 1] * b;
  }

  double T2a(0.0);
  for (int i = 0; i < N2a; ++i) {
    T2a += n[i] * apow[k[i]] * bpow[l[i]];
  }
  return T2a;
}


/* ******************************************************************
* Backward map for subregion 2b: formula 23.
****************************************************************** */
double
IAPWS97::Region2_SubregionFormulaB(double p, double h)
{
  const auto& k = inv2k_b;
  const auto& l = inv2l_b;
  const auto& n = inv2n_b;

  double pr = p;
  double hr = h / 2000.0;

  double a = pr - 2.0;
  double b = hr - 2.6;

  double apow[N2b_kmax];
  double bpow[N2b_lmax];

  apow[0] = 1.0;
  for (int i = 1; i < N2b_kmax; ++i) {
    apow[i] = apow[i - 1] * a;
  }

  bpow[0] = 1.0;
  for (int i = 1; i < N2b_lmax; ++i) {
    bpow[i] = bpow[i - 1] * b;
  }

  double T2b(0.0);
  for (int i = 0; i < N2b; ++i) {
    T2b += n[i] * apow[k[i]] * bpow[l[i]];
  }
  return T2b;
}


/* ******************************************************************
* Backward map for subregion 2c: formula 24.
****************************************************************** */
double
IAPWS97::Region2_SubregionFormulaC(double p, double h)
{
  const auto& k = inv2k_c;
  const auto& l = inv2l_c;
  const auto& n = inv2n_c;

  double pr = p;
  double hr = h / 2000.0;

  double a = pr + 25.0;
  double inva = 1.0 / a;
  double b = hr - 1.8;

  double apow[N2c_kmax];
  double bpow[N2c_lmax];
  double apow_inv[N2c_kmax_inv];

  apow[0] = 1.0;
  for (int i = 1; i < N2c_kmax; ++i) {
    apow[i] = apow[i - 1] * a;
  }

  apow_inv[0] = 1.0;
  for (int i = 1; i < N2c_kmax_inv; ++i) {
    apow_inv[i] = apow_inv[i - 1] * inva;
  }

  bpow[0] = 1.0;
  for (int i = 1; i < N2c_lmax; ++i) {
    bpow[i] = bpow[i - 1] * b;
  }

  double T2c(0.0), tmp;
  for (int i = 0; i < N2c; ++i) {
    tmp = (k[i] >= 0) ? apow[k[i]] : apow_inv[-k[i]];
    T2c += n[i] * tmp * bpow[l[i]];
  }
  return T2c;
}


/* ******************************************************************
* Region 3, Section 7, formula 28
****************************************************************** */
Properties
IAPWS97::Region3(double rho, double T)
{
  Properties prop;

  double rhor = rho / RHOC;
  double Tr = TC / T;

  double a = rhor;
  double b = Tr;

  double apow[N3_kmax];
  double bpow[N3_lmax];

  apow[0] = 1.0;
  for (int i = 1; i < N3_kmax; ++i)
    apow[i] = apow[i - 1] * a;

  bpow[0] = 1.0;
  for (int i = 1; i < N3_lmax; ++i)
    bpow[i] = bpow[i - 1] * b;

  const auto& k = rgn3k;
  const auto& l = rgn3l;
  const auto& n = rgn3n;

  double g(0.0), gd, gdd, gt(0.0), gtt(0.0), gdt(0.0), tmp;
  g = n[0] * std::log(a);
  gd = n[0] / a;
  gdd = -gd / a;

  for (int i = 1; i < N3; ++i) {
    tmp = bpow[l[i]];
    g += n[i] * apow[k[i]] * tmp;
    if (k[i] > 0) gd += n[i] * k[i] * apow[k[i] - 1] * tmp;
    if (k[i] > 1) gdd += n[i] * k[i] * (k[i] - 1) * apow[k[i] - 2] * tmp;

    tmp = apow[k[i]];
    if (l[i] > 0) gt += n[i] * l[i] * tmp * bpow[l[i] - 1];
    if (l[i] > 1) gtt += n[i] * l[i] * (l[i] - 1) * tmp * bpow[l[i] - 2];

    if (k[i] > 0 && l[i] > 0) gdt += n[i] * k[i] * l[i] * apow[k[i] - 1] * bpow[l[i] - 1];
  }

  double rhor2 = rhor * rhor;
  double A = rhor2 * (gd - Tr * gdt) * (gd - Tr * gdt);

  double p = rhor * gd * R * T * rho / 1000;
  prop.p = p;
  prop.T = T;
  prop.rho = rho;
  prop.v = 1.0 / rho;
  prop.h = R * T * (Tr * gt + rhor * gd);
  prop.u = R * T * Tr * gt; 
  prop.s = R * (Tr * gt - g);
  prop.cp = R * (-Tr * Tr * gtt + A / (2 * rhor * gd + rhor2 * gdd));
  prop.cv = -R * Tr * Tr * gtt;
  prop.w = std::sqrt(R * T * 1000 * (2 * rhor * gd + rhor2 * gdd - A / Tr / Tr / gtt));
  prop.kt = 1 / (2 * rhor * gd + rhor * rhor * gdd) / rho / R / T * 1000;
  prop.av = (gd - Tr * gdt) / (2 * gd + rhor * gdd) / T;
  prop.ap = (1 - Tr * gdt / gd) / T;
  prop.bp = rho * (2 + rhor * gdd / gd);

  prop.rgn = 3;
  prop.x = 1.0;
  if (T < TC && p < PC) {
    double Tsat = SaturationLineP(p);
    if (T < Tsat) prop.x = 0.0;
  }

  // derived properties
  prop.rho = 1.0 / prop.v;

  return prop;
}


/* ******************************************************************
* Backward map for Region 3: (p, T) -> specific volume
****************************************************************** */
double
IAPWS97::Region3_BackwardPT(double p, double T)
{
  auto id = Region3_SubregionIdPT(p, T);
  return Region3_SubregionFormula(p, T, id);
}


/* ******************************************************************
* Backward map for Region 3: (p, x)_T -> specific volume
****************************************************************** */
double
IAPWS97::Region3_BackwardPX(double p, double x, double T)
{
  auto id = Region3_SubregionIdPX(p, x);
  return Region3_SubregionFormula(p, T, id);
}


/* ******************************************************************
* Specific volume fitting
* Paper: "Revised Supplementary Release on Backward Equations for 
* Specific Volume as a Function of Pressure and Temperature v(p,T) 
* for Region 3 of the IAPWS Industrial Formulation 1997 for the 
* Thermodynamic Properties of Water and Steam
* http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Tables 2 and 10
****************************************************************** */
IAPWS97::ModelId
IAPWS97::Region3_SubregionIdPT(double p, double T)
{
  ModelId rgn;
  // top part of region 3 is plit into two parts
  if (p > 40.0) {
    rgn = (T <= Boundary3Formula2(p, rgn3k_ab, rgn3n_ab)) ? ModelId::a : ModelId::b;
  }

  // next part is plit into four parts along T-axis
  else if (25.0 < p  && p <= 40.0) {
    double Tcd = Boundary3Formula1(p, rgn3n_cd);
    double Tab = Boundary3Formula2(p, rgn3k_ab, rgn3n_ab);
    double Tef = Boundary3Formula3(p);

    if (T <= Tcd) {
      rgn = ModelId::c;
    } else if (Tcd < T && T <= Tab) {
      rgn = ModelId::d;
    } else if (Tab < T && T <= Tef) {
      rgn = ModelId::e;
    } else {
      rgn = ModelId::f;
    }
  }

  // next part is plit into six parts: c continued, d -> g and h, f -> i, j and k
  else if (23.5 < p && p <= 25.0) {
    double Tcd = Boundary3Formula1(p, rgn3n_cd);
    double Tgh = Boundary3Formula1(p, rgn3n_gh);
    double Tef = Boundary3Formula3(p);
    double Tij = Boundary3Formula1(p, rgn3n_ij);
    double Tjk = Boundary3Formula1(p, rgn3n_jk);

    if (T <= Tcd) {
      rgn = ModelId::c;
    } else if (Tcd < T && T <= Tgh) {
      rgn = ModelId::g;
    } else if (Tgh < T && T <= Tef) {
      rgn = ModelId::h;
    } else if (Tef < T && T <= Tij) {
      rgn = ModelId::i;
    } else if (Tij < T && T <= Tjk) {
      rgn = ModelId::j;
    } else {
      rgn = ModelId::k;
    }
  }

  // next part is plit into six parts: c, h, i, j, k are continued, g is continued as l
  else if (23.0 < p && p <= 23.5) {
    double Tcd = Boundary3Formula1(p, rgn3n_cd);
    double Tgh = Boundary3Formula1(p, rgn3n_gh);
    double Tef = Boundary3Formula3(p);
    double Tij = Boundary3Formula1(p, rgn3n_ij);
    double Tjk = Boundary3Formula1(p, rgn3n_jk);

    if (T <= Tcd) {
      rgn = ModelId::c;
    } else if (Tcd < T && T <= Tgh) {
      rgn = ModelId::l;
    } else if (Tgh < T && T <= Tef) {
      rgn = ModelId::h;
    } else if (Tef < T && T <= Tij) {
      rgn = ModelId::i;
    } else if (Tij < T && T <= Tjk) {
      rgn = ModelId::j;
    } else {
      rgn = ModelId::k;
    }
  }  

  // next part is split into eight parts: c, l, g, k are continues, h -> m and n, i -> o and p
  else if (22.5 < p && p <= 23.0) {
    double Tcd = Boundary3Formula1(p, rgn3n_cd);
    double Tgh = Boundary3Formula1(p, rgn3n_gh);
    double Tmn = Boundary3Formula1(p, rgn3n_mn);
    double Tef = Boundary3Formula3(p);
    double Top = Boundary3Formula2(p, rgn3k_op, rgn3n_op);
    double Tij = Boundary3Formula1(p, rgn3n_ij);
    double Tjk = Boundary3Formula1(p, rgn3n_jk);

    if (T <= Tcd) {
      rgn = ModelId::c;
    } else if (Tcd < T && T <= Tgh) {
      rgn = ModelId::l;
    } else if (Tgh < T && T <= Tmn) {
      rgn = ModelId::m;
    } else if (Tmn < T && T <= Tef) {
      rgn = ModelId::n;
    } else if (Tef < T && T <= Top) {
      rgn = ModelId::o;
    } else if (Top < T && T <= Tij) {
      rgn = ModelId::p;
    } else if (Tij < T && T <= Tjk) {
      rgn = ModelId::j;
    } else {
      rgn = ModelId::k;
    }
  }
  
  // next part is split into 10 parts: c and k are continued, rest are split differently
  else if (SaturationLineT(643.15) < p && p <= 22.5) {
    double Tcd = Boundary3Formula1(p, rgn3n_cd);
    double Tqu = Boundary3Formula1(p, rgn3n_qu);
    double Trx = Boundary3Formula1(p, rgn3n_rx);
    double Tjk = Boundary3Formula1(p, rgn3n_jk);

    if (T <= Tcd) {
      rgn = ModelId::c;
    } else if (Tcd < T && T <= Tqu) {
      rgn = ModelId::q;
    } else if (Tqu < T && T <= Trx) {
      // ranges from Table 10
      double Tef = Boundary3Formula3(p);
      double Twx = Boundary3Formula2(p, rgn3k_wx, rgn3n_wx);
      double Tuv = Boundary3Formula1(p, rgn3n_uv);
      if (22.11 < p && p <= 22.5) {
        if (T <= Tuv) {
          rgn = ModelId::u;
        } else if (Tuv <= T && T <= Tef) {
          rgn = ModelId::v;
        } else if (Tef <= T && T <= Twx) {
          rgn = ModelId::w;
        } else {
          rgn = ModelId::x;
        }
      } else if (22.064 < p && p <= 22.11) { // up to critical p
        if (T <= Tuv) {
          rgn = ModelId::u;
        } else if (Tuv <= T && T <= Tef) {
          rgn = ModelId::y;
        } else if (Tef <= T && T <= Twx) {
          rgn = ModelId::z;
        } else {
          rgn = ModelId::x;
        }
      } else if (T > SaturationLineP(p)) {
        if (SaturationLineT(643.15) < p && p <= 21.90096265) {
          rgn = ModelId::x;
        } else if (21.90096265 < p && p <= 22.064) {
          if (T <= Twx) rgn = ModelId::z;
          else rgn = ModelId::x;
        }
      } else if (T <= SaturationLineP(p)) {
        if (SaturationLineT(643.15) < p && p <= 21.93161551) {
          rgn = ModelId::u;
        } else if (21.93161551 < p && p <= 22.064) {
          if (T <= Tuv) rgn = ModelId::u;
          else rgn = ModelId::y;
        }
      }
    } else if (Trx < T && T <= Tjk) {
      rgn = ModelId::r;
    } else {
      rgn = ModelId::k;
    }
  }

  // next part is split into 4 parts: c, r and k are continued; q is continued as s
  else if (20.5 < p && p <= SaturationLineT(643.15)) {
    double Tcd = Boundary3Formula1(p, rgn3n_cd);
    double Ts = SaturationLineP(p);
    double Tjk = Boundary3Formula1(p, rgn3n_jk);
    if (T <= Tcd) {
      rgn = ModelId::c;
    } else if (Tcd < T && T <= Ts) {
      rgn = ModelId::s;
    } else if (Ts < T && T <= Tjk) {
      rgn = ModelId::r;
    } else {
      rgn = ModelId::k;
    }
  }

  // next part is split into 3 parts: c and s are continued; k is contiuned as t
  else if (19.00881189173929 < p && p <= 20.5) {
    double Tcd = Boundary3Formula1(p, rgn3n_cd);
    double Ts = SaturationLineP(p);
    if (T <= Tcd) {
      rgn = ModelId::c;
    } else if (Tcd < T && T <= Ts) {
      rgn = ModelId::s;
    } else {
      rgn = ModelId::t;
    }
  }

  // final part is split into 2 parts: c and t are continued
  else if (PSAT_623 < p && p <= 19.00881189173929) {
    double Ts = SaturationLineP(p);
    if (T <= Ts) rgn = ModelId::c;
    else rgn = ModelId::t;
  }

  return rgn;
}


/* ******************************************************************
*
****************************************************************** */
IAPWS97::ModelId
IAPWS97::Region3_SubregionIdPX(double p, double x)
{
  ModelId rgn;

  // saturated liquid
  if (x == 0.0) {
    if (p < 19.00881189) {
      rgn = ModelId::c;
    } else if (p < 21.0434) {
      rgn = ModelId::s;
    } else if (p < 21.9316) {
      rgn = ModelId::u;
    } else {
      rgn = ModelId::y;
    }
  }
  // saturated vapor
  else {
    if (p < 20.5) {
      rgn = ModelId::t;
    } else if (p < 21.0434) {
      rgn = ModelId::r;
    } else if (p < 21.9009) {
      rgn = ModelId::x;
    } else {
      rgn = ModelId::z;
    }
  }

  return rgn;
}


/* ******************************************************************
* Backward map for Region 3 return specific volume, m3/kg
* http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, formulas 4 and 5
****************************************************************** */
double
IAPWS97::Region3_SubregionFormula(double p, double T, IAPWS97::ModelId id)
{
  const auto& d = inv3_table_d[static_cast<int>(id)];
  const auto& k = d.k;
  const auto& l = d.l;
  const auto& n = d.n;
  const auto& par = d.par;

  double pr = p / par[1];
  double Tr = T / par[2];

  double a = pr - par[3];
  double b = Tr - par[4];
  if (id != IAPWS97::ModelId::n) {
    a = std::pow(a, par[5]);
    b = std::pow(b, par[6]);
  }

  double inva = 1.0 / a;
  double invb = 1.0 / b;

  // compute extrema
  int size = d.size;
  int kneg(0), kpos(0), lneg(0), lpos(0);
  for (int i = 0; i < size; ++i) {
    kneg = std::min(kneg, k[i]);
    kpos = std::max(kpos, k[i]);

    lneg = std::min(lneg, l[i]);
    lpos = std::max(lpos, l[i]);
  }

  double apow[N3_max], apow_inv[N3_max];
  double bpow[N3_max], bpow_inv[N3_max];

  apow[0] = 1.0;
  for (int i = 1; i <= kpos; ++i)
    apow[i] = apow[i - 1] * a;

  apow_inv[0] = 1.0;
  for (int i = 1; i <= -kneg; ++i)
    apow_inv[i] = apow_inv[i - 1] * inva;

  bpow[0] = 1.0;
  for (int i = 1; i <= lpos; ++i)
    bpow[i] = bpow[i - 1] * b;

  bpow_inv[0] = 1.0;
  for (int i = 1; i <= -lneg; ++i)
    bpow_inv[i] = bpow_inv[i - 1] * invb;

  double v(0.0), tmpk, tmpl;
  for (int i = 0; i < size; ++i) {
    tmpk = (k[i] >= 0) ? apow[k[i]] : apow_inv[-k[i]];
    tmpl = (l[i] >= 0) ? bpow[l[i]] : bpow_inv[-l[i]];
    v += n[i] * tmpk * tmpl;
  }

  if (id == IAPWS97::ModelId::n) v = std::exp(v);
  else v = std::pow(v, par[7]);

  return par[0] * v;
}


/* ******************************************************************
* Backward map for Region 3: (p, h) -> v
* doi: 10.1007/978-3-540-74234-0, formulas (2.26)-(2.27)
****************************************************************** */
double
IAPWS97::Region3_BackwardPH_v(double p, double h)
{
  double hf = BoundaryLine3a3b(p);
  if (h <= hf) {
    return Region3_SubregionFormulaA_v(p, h);
  }
  return Region3_SubregionFormulaB_v(p, h);
}


/* ******************************************************************
* Backward map for subregion 3a: formulas (2.26)
****************************************************************** */
double
IAPWS97::Region3_SubregionFormulaA_v(double p, double h)
{
  const auto& k = inv3ak_ph_v;
  const auto& l = inv3al_ph_v;
  const auto& n = inv3an_ph_v;

  double pr = p / 100.0;
  double hr = h / 2100.0;

  double a = pr + 0.128;
  double b = hr - 0.727;
  double inva = 1.0 / a;

  double apow[N3av_kmax];
  double bpow[N3av_lmax];
  double apow_inv[N3av_kmax_inv];

  apow[0] = 1.0;
  for (int i = 1; i < N3av_kmax; ++i) {
    apow[i] = apow[i - 1] * a;
  }

  apow_inv[0] = 1.0;
  for (int i = 1; i < N3av_kmax_inv; ++i) {
    apow_inv[i] = apow_inv[i - 1] * inva;
  }

  bpow[0] = 1.0;
  for (int i = 1; i < N3av_lmax; ++i) {
    bpow[i] = bpow[i - 1] * b;
  }

  double v(0.0), tmp;
  for (int i = 0; i < N3av; ++i) {
    tmp = (k[i] >= 0) ? apow[k[i]] : apow_inv[-k[i]];
    v += n[i] * tmp * bpow[l[i]];
  }
  return 0.0028 * v;
}


/* ******************************************************************
* Backward map for subregion 3b: formulas (2.27)
****************************************************************** */
double
IAPWS97::Region3_SubregionFormulaB_v(double p, double h)
{
  const auto& k = inv3bk_ph_v;
  const auto& l = inv3bl_ph_v;
  const auto& n = inv3bn_ph_v;

  double pr = p / 100.0;
  double hr = h / 2800.0;

  double a = pr + 0.0661;
  double b = hr - 0.720;
  double inva = 1.0 / a;

  double apow[N3bv_kmax];
  double bpow[N3bv_lmax];
  double apow_inv[N3bv_kmax_inv];

  apow[0] = 1.0;
  for (int i = 1; i < N3bv_kmax; ++i) {
    apow[i] = apow[i - 1] * a;
  }

  apow_inv[0] = 1.0;
  for (int i = 1; i < N3bv_kmax_inv; ++i) {
    apow_inv[i] = apow_inv[i - 1] * inva;
  }

  bpow[0] = 1.0;
  for (int i = 1; i < N3bv_lmax; ++i) {
    bpow[i] = bpow[i - 1] * b;
  }

  double v(0.0), tmp;
  for (int i = 0; i < N3bv; ++i) {
    tmp = (k[i] >= 0) ? apow[k[i]] : apow_inv[-k[i]];
    v += n[i] * tmp * bpow[l[i]];
  }
  return 0.0088 * v;
}


/* ******************************************************************
* Backward map for Region 3: (p, h) -> T
* doi: 10.1007/978-3-540-74234-0, formulas (2.28)-(2.29)
****************************************************************** */
double
IAPWS97::Region3_BackwardPH_T(double p, double h)
{
  double hf = BoundaryLine3a3b(p);
  if (h <= hf) {
    return Region3_SubregionFormulaA_T(p, h);
  }
  return Region3_SubregionFormulaB_T(p, h);
}


/* ******************************************************************
* Backward map for subregion 3a: formulas (2.28)
****************************************************************** */
double
IAPWS97::Region3_SubregionFormulaA_T(double p, double h)
{
  const auto& k = inv3ak_ph_T;
  const auto& l = inv3al_ph_T;
  const auto& n = inv3an_ph_T;

  double pr = p / 100.0;
  double hr = h / 2300.0;

  double a = pr + 0.240;
  double b = hr - 0.615;
  double inva = 1.0 / a;

  double apow[N3aT_kmax];
  double bpow[N3aT_lmax];
  double apow_inv[N3aT_kmax_inv];

  apow[0] = 1.0;
  for (int i = 1; i < N3aT_kmax; ++i) {
    apow[i] = apow[i - 1] * a;
  }

  apow_inv[0] = 1.0;
  for (int i = 1; i < N3aT_kmax_inv; ++i) {
    apow_inv[i] = apow_inv[i - 1] * inva;
  }

  bpow[0] = 1.0;
  for (int i = 1; i < N3aT_lmax; ++i) {
    bpow[i] = bpow[i - 1] * b;
  }

  double T(0.0), tmp;
  for (int i = 0; i < N3aT; ++i) {
    tmp = (k[i] >= 0) ? apow[k[i]] : apow_inv[-k[i]];
    T += n[i] * tmp * bpow[l[i]];
  }
  return 760.0 * T;
}


/* ******************************************************************
* Backward map for subregion 3b: formulas (2.29)
****************************************************************** */
double
IAPWS97::Region3_SubregionFormulaB_T(double p, double h)
{
  const auto& k = inv3bk_ph_T;
  const auto& l = inv3bl_ph_T;
  const auto& n = inv3bn_ph_T;

  double pr = p / 100.0;
  double hr = h / 2800.0;

  double a = pr + 0.298;
  double b = hr - 0.720;
  double inva = 1.0 / a;

  double apow[N3bT_kmax];
  double bpow[N3bT_lmax];
  double apow_inv[N3bT_kmax_inv];

  apow[0] = 1.0;
  for (int i = 1; i < N3bT_kmax; ++i) {
    apow[i] = apow[i - 1] * a;
  }

  apow_inv[0] = 1.0;
  for (int i = 1; i < N3bT_kmax_inv; ++i) {
    apow_inv[i] = apow_inv[i - 1] * inva;
  }

  bpow[0] = 1.0;
  for (int i = 1; i < N3bT_lmax; ++i) {
    bpow[i] = bpow[i - 1] * b;
  }

  double T(0.0), tmp;
  for (int i = 0; i < N3bT; ++i) {
    tmp = (k[i] >= 0) ? apow[k[i]] : apow_inv[-k[i]];
    T += n[i] * tmp * bpow[l[i]];
  }
  return 860.0 * T;
}


/* ******************************************************************
* Region 4, (p, x) 
****************************************************************** */
Properties
IAPWS97::Region4(double p, double x)
{
  double T, rhol, rhov;
  Properties prop, prop1, prop2;

  T = SaturationLineP(p);
  if (T > 623.15) {
    rhol = 1.0 / Region3_BackwardPX(p, 0.0, T);
    prop1 = Region3(rhol, T);

    rhov = 1.0 / Region3_BackwardPX(p, 1.0, T);
    prop2 = Region3(rhov, T);
  } else {
    prop1 = Region1(p, T);
    prop2 = Region2(p, T);
  }

  prop.T = T;
  prop.p = p;
  prop.v = (1.0 - x) * prop1.v + x * prop2.v;
  prop.h = (1.0 - x) * prop1.h + x * prop2.h;
  prop.s = (1.0 - x) * prop1.s + x * prop2.s;
  prop.x = x;
  prop.rgn = 4;

  // derived properties
  prop.rho = 1.0 / prop.v;

  return prop;
}


/* ******************************************************************
* Region 5, Section 9, formulas 33 and 34
****************************************************************** */
Properties
IAPWS97::Region5(double p, double T)
{
  Properties prop;

  double pr = p;
  double Tr = 1000.0 / T;
  double a = pr;
  double b = Tr;

  // ideal gas part of the dimensionless Gibbs free energy
  const auto& no = rgn5n_o;

  double go, gop, gopp, got, gott, gopt(0.0);
  go = std::log(pr);
  gop = 1.0 / pr;
  gopp = -gop / pr;
  
  double b2, bneg1, bneg2, bneg3, bneg4, bneg5;
  b2 = b * b;
  bneg1 = 1.0 / b;
  bneg2 = bneg1 * bneg1;
  bneg3 = bneg2 * bneg1;
  bneg4 = bneg3 * bneg1;
  bneg5 = bneg4 * bneg1;

  go += no[0] + no[1] * b + no[2] * bneg3 + no[3] * bneg2 + no[4] * bneg1 + no[5] * b2;

  got = no[1] - 3 * no[2] * bneg4 - 2 * no[3] * bneg3 - no[4] * bneg2 + 2 * no[5] * b;

  gott = 12 * no[2] * bneg5 + 6 * no[3] * bneg4 + 2 * no[4] * bneg3 + 2 * no[5];

  // residual part of the Gibbs free energy
  double apow[N5r_kmax];
  double bpow[N5r_lmax];

  apow[0] = 1.0;
  for (int i = 1; i < N5r_kmax; ++i)
    apow[i] = apow[i - 1] * a;

  bpow[0] = 1.0;
  for (int i = 1; i < N5r_lmax; ++i)
    bpow[i] = bpow[i - 1] * b;
 
  const auto& k = rgn5k_r;
  const auto& l = rgn5l_r;
  const auto& nr = rgn5n_r;

  double gr(0.0), grp(0.0), grpp(0.0), grt(0.0), grtt(0.0), grpt(0.0), tmp;
  for (int i = 0; i < N5r; ++i) {
    tmp = bpow[l[i]];
    gr += nr[i] * apow[k[i]] * tmp;
    if (k[i] > 0) grp += nr[i] * k[i] * apow[k[i] - 1] * tmp;
    if (k[i] > 1) grpp += nr[i] * k[i] * (k[i] - 1) * apow[k[i] - 2] * tmp;

    tmp = apow[k[i]];
    if (l[i] > 0) grt += nr[i] * l[i] * tmp * bpow[l[i] - 1];
    if (l[i] > 1) grtt += nr[i] * l[i] * (l[i] - 1) * tmp * bpow[l[i] - 2];

    if (k[i] > 0 && l[i] > 0) grpt += nr[i] * k[i] * l[i] * apow[k[i] - 1] * bpow[l[i] - 1];
  }

  double A = (1 + pr * grp - Tr * pr * grpt) * (1 + pr * grp - Tr * pr * grpt);
  double B = 1 + 2 * pr * grp + pr * pr * grp * grp;
  double C = ((gop + grp) - Tr * (gopt + grpt)) * ((gop + grp) - Tr * (gopt + grpt));

  prop.p = p;
  prop.T = T;
  prop.v = (pr / p) * (gop + grp) * R * T / 1000;
  prop.h = Tr * (got + grt) * R * T;
  prop.u = R * T * (Tr * (got + grt) - pr * (gop + grp));
  prop.s = R * (Tr * (got + grt) - (go + gr));
  prop.cp = -R * Tr * Tr * (gott + grtt);
  // prop.cv = R * (-Tr * Tr * (gott + grtt) + C / (gopp + grpp)); (similar to below)
  prop.cv = R * (-Tr * Tr * (gott + grtt) - A / (1 - pr * pr * grpp));
  prop.w = std::sqrt(R * T * 1000 * B / (1 - pr * pr * grpp + A / Tr / Tr / (gott + grtt)));
  prop.kt = (1 - pr * pr * grpp) / (1 + pr * grp) / p;
  prop.av = (1 + pr * grp - Tr * pr * grpt) / (1 + pr * grp) / T;
  prop.x = 1.0;
  prop.rgn = 5;

  // derived properties
  prop.rho = 1.0 / prop.v;

  return prop;
}


/* ******************************************************************
* Region4, Section 8.2, formula 31
* Validity range: 611.213 Pa < p < 22.064 MPa
****************************************************************** */
double
IAPWS97::SaturationLineP(double p)
{
  const auto& n1 = rgn4n1;
  const auto& n2 = rgn4n2;
  const auto& n3 = rgn4n3;
  const auto& n4 = rgn4n4;
  const auto& n5 = rgn4n5;
  const auto& n6 = rgn4n6;
  const auto& n7 = rgn4n7;
  const auto& n8 = rgn4n8;
  const auto& n9 = rgn4n9;
  const auto& n10 = rgn4n10;

  double beta, beta2, E, F, G, D;
  beta = std::pow(p, 0.25);
  beta2 = beta * beta;

  E = beta2 + n3 * beta + n6;
  F = n1 * beta2 + n4 * beta + n7;
  G = n2 * beta2 + n5 * beta + n8;

  D = 2 * G / (-F - std::pow((F * F - 4 * E * G), 0.5));
  return (n10 + D - std::pow(((n10 + D) * (n10 + D) - 4 * (n9 + n10 * D)), 0.5)) / 2;
}


/* ******************************************************************
* Region4, Section 8.2, formula 30
* Validity range: 273.15 K <= T <= 647.096 K
****************************************************************** */
double
IAPWS97::SaturationLineT(double T)
{
  const auto& n1 = rgn4n1;
  const auto& n2 = rgn4n2;
  const auto& n3 = rgn4n3;
  const auto& n4 = rgn4n4;
  const auto& n5 = rgn4n5;
  const auto& n6 = rgn4n6;
  const auto& n7 = rgn4n7;
  const auto& n8 = rgn4n8;
  const auto& n9 = rgn4n9;
  const auto& n10 = rgn4n10;

  double tau, tau2, A, B, C;
  tau = T + n9 / (T - n10);
  tau2 = tau * tau;

  A = tau2  + n1 * tau + n2;
  B = n3 * tau2 + n4 * tau + n5;
  C = n6 * tau2 + n7 * tau + n8;
  return std::pow(2 * C / (-B + std::sqrt(B * B - 4 * A * C)), 4);
}


/* ******************************************************************
* Boundary between regions 2 and 3, Section 4, formula 6
****************************************************************** */
double
IAPWS97::BoundaryLine23(double p)
{
  const auto& n3 = bnd23_n3;
  const auto& n4 = bnd23_n4;
  const auto& n5 = bnd23_n5;
  return n4 + std::pow((p - n5) / n3, 0.5);
}


/* ******************************************************************
* Boundary between regions 3e and 3f, formula 3
****************************************************************** */
double
IAPWS97::Boundary3Formula3(double p)
{
  return 3.727888004 * (p - 22.064) + 647.096;
}


/* ******************************************************************
* Boundary between regions 2b and 2c, formula 21
****************************************************************** */
double
IAPWS97::BoundaryLine2b2c(double p)
{
  return 0.26526571908428e4 + std::sqrt((p - 4.5257578905948) / 1.2809002730136e-4);
}


/* ******************************************************************
* Boundary between regions 3a and 3b,
* doi: 10.1007/978-3-540-74234-0, formula (2.25)
****************************************************************** */
double
IAPWS97::BoundaryLine3a3b(double p)
{
  double h = 0.201464004206875e4 + 3.74696550136983 * p 
           - 0.0219921901054187 * p * p 
           + 0.875131686009950e-4 * p * p * p;
  return h;
}


/* ******************************************************************
* Boundary of region 3,
* doi: 10.1007/978-3-540-74234-0, formula (2.18)
****************************************************************** */
double
IAPWS97::SaturationLineH(double h)
{
  const auto& k = rgn4k_h;
  const auto& l = rgn4l_h;
  const auto& n = rgn4n_h;

  double hr = h / 2600.0;
  double a = hr - 1.02;
  double b = hr - 0.608;

  double apow[N4_kmax];
  double bpow[N4_lmax];

  apow[0] = 1.0;
  for (int i = 1; i < N4_kmax; ++i) {
    apow[i] = apow[i - 1] * a;
  }

  bpow[0] = 1.0;
  for (int i = 1; i < N4_lmax; ++i) {
    bpow[i] = bpow[i - 1] * b;
  }

  double tmp(0.0);
  for (int i = 0; i < N4; ++i) {
    tmp += n[i] * apow[k[i]] * bpow[l[i]];
  }
  return 22 * tmp;
}


/* ******************************************************************
* Map P/T to region Id
****************************************************************** */
int 
IAPWS97::RegionIdPT(double p, double T)
{
  int rgn(0);

  if (1073.15 < T && T <= 2273.15 && PMIN <= p && p <= 50.0) {
    rgn = 5;
  } else if (PMIN <= p && p <= PSAT_623) {
    double Tsat = SaturationLineP(p);
    if (273.15 <= T && T <= Tsat) {
      rgn = 1;
    } else if (Tsat < T && T <= 1073.15) {
      rgn = 2;
    }
  } else if (PSAT_623 < p && p <= 100.0) {
    double Tbnd = BoundaryLine23(p);
    if (273.15 <= T && T <= 623.15) {
      rgn = 1;
    } else if (623.15 < T && T < Tbnd) {
      rgn = 3;
    } else if (Tbnd <= T && T <= 1073.15) {
      rgn = 2;
    }
  }
  return rgn;
}


/* ******************************************************************
* Map P/H to region Id
* doi: 10.1007/978-3-540-74234-0. Fig. 2.5
****************************************************************** */
int
IAPWS97::RegionIdPH(double p, double h)
{
  int rgn(0);

  double hmin = (Region1(p, 273.15)).h;
  double hmax = (Region5(p, 2273.15)).h;

  if (PMIN <= p && p <= PSAT_623) {
    double Tsat = SaturationLineP(p);
    double h14 = (Region1(p, Tsat)).h;
    double h24 = (Region2(p, Tsat)).h;
    double h25 = (Region2(p, 1073.15)).h;
    if (hmin <= h && h <= h14) {
      rgn = 1;
    } else if (h14 < h && h < h24) {
      rgn = 4;
    } else if (h24 <= h && h <= h25) {
      rgn = 2;
    } else if (h25 < h && h <= hmax) {
      rgn = 5;
    }
  } else if (PSAT_623 < p && p < PC) {
    double h13 = (Region1(p, 623.15)).h;
    double h32 = (Region2(p, BoundaryLine23(p))).h;
    double h25 = (Region2(p, 1073.15)).h;
    if (hmin <= h && h <= h13) {
      rgn = 1;
    } else if (h13 < h && h < h32) {
      double p34;
      double hmin3 = (Region1(PSAT_623, 623.15)).h;
      double hmax3 = (Region2(PSAT_623, 623.15)).h;
      if (hmin3 <= h && h <= hmax3) {
        p34 = SaturationLineH(h);
      } else {
        p34 = PSAT_623;
      }
      if (p < p34) {
        rgn = 4;
      } else {
        rgn = 3;
      }
    } else if (h32 <= h && h <= h25) {
      rgn = 2;
    } else if (h25 < h && h <= hmax) {
      rgn = 5;
    }
  } else if (PC <= p && p <= 100.0) {
    double h13 = (Region1(p, 623.15)).h;
    double h32 = (Region2(p, BoundaryLine23(p))).h;
    double h25 = (Region2(p, 1073.15)).h;
    if (hmin <= h && h <= h13) {
      rgn = 1;
    } else if (h13 < h && h < h32) {
      rgn = 3;
    } else if (h32 <= h && h <= h25) {
      rgn = 2;
    } else if (p <= 50.0 && h25 <= h && h <= hmax) {
      rgn = 5;
    }
  }

  return rgn;
}


/* ******************************************************************
* Phase
****************************************************************** */
Phase_t
IAPWS97::PhaseId(double p, double T, int rgn, double x)
{
  Phase_t phase;

  if (p > PC && T > TC) {
    phase = Phase_t::SupercriticalLiquid;
  } else if (T > TC) {
    phase = Phase_t::Gas;
  } else if (p > PC) {
    phase = Phase_t::CompressibleLiquid;
  } else if (p == PC && T == TC) {
    phase = Phase_t::CriticalPoint;
  } else if (rgn == 4 && x == 1.0) {
    phase = Phase_t::SaturatedVapor;
  } else if (rgn == 4 && x == 0.0) {
    phase = Phase_t::SaturatedLiquid;
  } else if (rgn == 4) {
    phase = Phase_t::TwoPhases;
  } else if (x == 1.0) {
    phase = Phase_t::Vapor;
  } else if (x == 0.0) {
    phase = Phase_t::Liquid;
  }

  return phase;
}


/* ******************************************************************
* Finalize properties
****************************************************************** */
Properties
IAPWS97::ExtendProperies(const Properties& prop_in)
{
  Properties prop;
  prop = prop_in;

  double p = prop.p;
  double T = prop.T;
  prop.helmholtz = prop.u - T * prop.s;
  prop.gibbs = prop.h - T * prop.s;

  prop.mu = Viscosity(prop.rho, T);
  prop.k = ThermalConductivity(prop.rho, T, prop);

  return prop;
}


/* ******************************************************************
* Thermal conductivity
* http://www.iapws.org/relguide/ThCond.html, formulas (15)-(17)
****************************************************************** */
double
IAPWS97::ThermalConductivity(double rho, double T, Properties& prop)
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

  // critical enhancement (k2) is not implemented yet FIXME
  return 1e-3 * k0 * k1;
}


/* ******************************************************************
* Surface tension
* http://www.iapws.org/relguide/Surf-H2O.html
* Validity range 248.15 <= T <= TC
****************************************************************** */
double
IAPWS97::SurfaceTension(double T)
{
  if (248.15 <= T && T <= TC) {
    double tau = 1.0 - T / TC;
    return 235.8e-3 * std::pow(tau, 1.256) * (1.0 - 0.625 * tau);
  }
  return -1.0;
}


/* ******************************************************************
* Dynamic viscosity
* http://www.iapws.org/relguide/viscosity.html, formulas (10)-(12)
****************************************************************** */
double
IAPWS97::Viscosity(double rho, double T)
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

  // third factor (critical enhancement) is not implemeneted yet FIXME

  return 1e-6 * mu0 * mu1;
}


/* ******************************************************************
* i/o
****************************************************************** */
void
IAPWS97::Print(Properties& prop)
{
  std::cout << std::setprecision(12)
    << "========  region: " << prop.rgn << "  =========" 
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
    << "\nappha_p = " << prop.ap
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

