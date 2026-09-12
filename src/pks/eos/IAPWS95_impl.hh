/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_IAPWS95_IMPL_HH_
#define AMANZI_IAPWS95_IMPL_HH_

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* Residual part of Hemholtz free energy
* http://www.iapws.org/relguide/IAPWS-95.html
****************************************************************** */
template<bool SecondOrder>
auto IAPWS95::ResidualPartImpl_(double rho, double T)
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

    if constexpr (SecondOrder) {
      if (t1[i] != 1.0) gtt += n1[i] * t1[i] * (t1[i] - 1.0) * tmp1 * std::pow(tau, t1[i] - 2.0);
      if (d1[i] > 0) gdt += n1[i] * d1[i] * t1[i] * dpow[d1[i] - 1] * tmp2;
    }
  }

  // exponential terms
  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
  for (int i = 0; i < 44; ++i) {
    int c = c2[i];
    int d = d2[i];
    int t = t2[i];
    double tmp1 = dpow[d] * epow[c];
    double tmp2 = tpow[t];

    g += n2[i] * tmp1 * tmp2;
    gd += n2[i] * tmp2 * dpow[d - 1] * (d - c * dpow[c]) * epow[c]; 
    gt += n2[i] * t * tmp1 * tpow[t - 1];

    if constexpr (SecondOrder) {
      if (d == 1) {
        gdd += n2[i] * tmp2 * epow[c] * c * dpow[c - 1] * (-1 + c * dpow[c] - c);
      } else {
        gdd += n2[i] * tmp2 * epow[c] * dpow[d - 2] * ((d - c * dpow[c]) * (d - 1 - c * dpow[c]) - c * c * dpow[c]);
      }

      if (t > 1) gtt += n2[i] * t * (t - 1) * tmp1 * tpow[t - 2];

      gdt += n2[i] * t * tpow[t - 1] * dpow[d - 1] * (d - c * dpow[c]) * epow[c]; 
    }
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

    g += n3[i] * tmp4 * tmp3;

    // we know that d > 1, so the result is "symmetric" to gt and gtt
    gd += n3[i] * tmp3 * dpow[d - 1] * tpow[t] * (d - 2 * al * delta * (delta - 1));
    gt += n3[i] * tmp4 * tmp3 * (t / tau - 2 * be * (tau - ga));

    if constexpr (SecondOrder) {
      tmp5 = t / tau - 2 * be * (tau - ga);

      gdd += n3[i] * dpow[d - 2] * tpow[t] * tmp3 *
        (2 * al * delta * delta * (2 * tmp1 - 1.0) - 4 * d * al * delta * (delta - 1) + d * (d - 1));

      gtt += n3[i] * tmp4 * tmp3 * (tmp5 * tmp5 - t / tau / tau - 2 * be);

      gdt += n3[i] * tmp4 * tmp3 * (d / delta - 2 * al * (delta - 1)) * (t / tau - 2 * be * (tau - ga)); 
    }
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

    if (delta != 1.0) {
      tmp5 = theta * A[i] / beta4[i] * tmp3 / tmp1 + B[i] * a4[i] * tmp4 / tmp1;
      tmp6 = A[i] / beta4[i] * tmp3 / tmp1;
      deld = 2 * (delta - 1) * tmp5;
    }

    delt = -2 * theta;

    tmp7 = std::pow(del, b4[i]);
    tmp8 = b4[i] * std::pow(del, b4[i] - 1);
    g += n4[i] * tmp7 * delta * psi;

    gd += n4[i] * (tmp7 * (psi + delta * psid) + tmp8 * deld * delta * psi);
    gt += n4[i] * delta * (psit * tmp7 + tmp8 * psi * delt);

    if constexpr (SecondOrder) {
      tmp9 = b4[i] * (b4[i] - 1) * std::pow(del, b4[i] - 2);

      if (delta != 1.0) {
        deldd = 2 * tmp5 + tmp1 * (4 * B[i] * a4[i] * (a4[i] - 1) * tmp4 / tmp1 / tmp1 +
                                   2 * tmp6 * tmp6 + 
                                   4 * A[i] * theta / beta4[i] * (0.5 / beta4[i] - 1) * tmp3 / tmp1 / tmp1);

        deldt = -2 * tmp6 * (delta - 1);
      }
      deltt = 2.0;

      psidd = 2 * C[i] * (2 * C[i] * tmp1 - 1) * psi;
      psitt = 2 * D[i] * (2 * D[i] * tmp2 - 1) * psi;
      psidt = 4 * C[i] * D[i] * (delta - 1) * (tau - 1) * psi;

      gdd += n4[i] * (tmp7 * (2 * psid + delta * psidd) + 
                      2 * tmp8 * deld * (psi + delta * psid) +
                      (tmp8 * deldd + tmp9 * deld * deld) * delta * psi);

      gtt += n4[i] * delta * (tmp7 * psitt + 2 * tmp8 * psit * delt 
                                           + tmp9 * delt * delt * psi + tmp8 * deltt * psi);
    
      gdt += n4[i] * (tmp7 * (psit + delta * psidt) + delta * tmp8 * deld * psit + 
                      tmp8 * delt * (psi + delta * psid) + 
                      (tmp9 * deld * delt + tmp8 * deldt) * delta * psi);
    }
  }

  residual_calls++;

  if constexpr (SecondOrder)
    return std::array<double, 6>{ g, gd, gt, gdd, gdt, gtt };
  else
    return std::array<double, 3>{ g, gd, gt };
}


/* ******************************************************************
* Specializations
****************************************************************** */
inline
std::array<double, 3>
IAPWS95::ResidualPartFirst(double rho, double T) {
  return ResidualPartImpl_<false>(rho, T);
}

inline
std::array<double, 6>
IAPWS95::ResidualPart(double rho, double T) {
  return ResidualPartImpl_<true>(rho, T);
}

} // namespace AmanziEOS
} // namespace Amanzi

#endif
