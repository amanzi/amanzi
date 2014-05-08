/*

Functor for evaluating QSat

Author: Ethan Coon (ecoon@lanl.gov)
        Chonggang Xu (cxu@lanl.gov)

Licencse: BSD
*/

#include <cmath>
#include <algorithm>
#include "vegetation.hh"

namespace Amanzi {
namespace BGC {
namespace BGCSimple {

QSat::QSat() :
    a0(6.11213476),
	a1(0.444007856),
	a2(0.143064234 * std::pow(10.0, -1)),
	a3(0.264461437 * std::pow(10.0, -3)),
	a4(0.305903558 * std::pow(10.0, -5)),
	a5(0.196237241 * std::pow(10.0, -7)),
	a6(0.892344772 * std::pow(10.0, -10)),
	a7(-0.373208410 * std::pow(10.0, -12)),
	a8(0.209339997 * std::pow(10.0, -15)),
	b0(0.444017302),
	b1(0.286064092 * std::pow(10.0, -1)),
	b2(0.794683137 * std::pow(10.0, -3)),
	b3(0.121211669 * std::pow(10.0, -4)),
	b4(0.103354611 * std::pow(10.0, -6)),
	b5(0.404125005 * std::pow(10.0, -9)),
	b6(-0.788037859 * std::pow(10.0, -12)),
	b7(-0.114596802 * std::pow(10.0, -13)),
	b8(0.381294516 * std::pow(10.0, -16)),
	c0(6.11123516),
	c1(0.503109514),
	c2(0.188369801 * std::pow(10.0, -1)),
	c3(0.420547422 * std::pow(10.0, -3)),
	c4(0.614396778 * std::pow(10.0, -5)),
	c5(0.602780717 * std::pow(10.0, -7)),
	c6(0.387940929 * std::pow(10.0, -9)),
	c7(0.149436277 * std::pow(10.0, -11)),
	c8(0.262655803 * std::pow(10.0, -14)),
	d0(0.503277922),
	d1(0.377289173 * std::pow(10.0, -1)),
	d2(0.126801703 * std::pow(10.0, -2)),
	d3(0.249468427 * std::pow(10.0, -4)),
	d4(0.313703411 * std::pow(10.0, -6)),
	d5(0.257180651 * std::pow(10.0, -8)),
	d6(0.133268878 * std::pow(10.0, -10)),
	d7(0.394116744 * std::pow(10.0, -13)),
	d8(0.498070196 * std::pow(10.0, -16))
{}


void
QSat::operator()(double tleafk, double pressure,
                 double* es, double* esdT, double* qs, double* qsdT) {
  double td = tleafk - 273.15;
  td = std::min(std::max(-75.0, td), 100.0);

  if (td >= 0.0) {
    *es = a0 + td * (a1 + td * (a2 + td * (a3 + td * (a4 + td * (a5 + td * (a6 + td * (a7 + td * a8)))))));
    *esdT = b0 + td * (b1 + td * (b2 + td * (b3 + td * (b4 + td * (b5 + td * (b6 + td * (b7 + td * b8)))))));
  } else {
    *es = c0 + td * (c1 + td * (c2 + td * (c3 + td * (c4 + td * (c5 + td * (c6 + td * (c7 + td * c8)))))));
    *esdT = d0 + td * (d1 + td * (d2 + td * (d3 + td * (d4 + td * (d5 + td * (d6 + td * (d7 + td * d8)))))));
  }

  *es = *es * 100; // [Pa]
  *esdT = *esdT * 100; // [Pa/K]

  double vp = 1.0 / (pressure - 0.378 * (*es));
  double vp1 = 0.622 * vp;
  double vp2 = vp1 * vp;
  *qs = (*es) * vp1;  // [kg/kg]
  *qsdT = (*esdT) * vp2 * pressure; // [1/K]
  return;
}


double DayLength(double lat, int doy) {
  const double PI = 3.141592653589793;
  double LatRad = lat * (2.0 * PI) / 360.0;
  LatRad = std::max(std::min(LatRad, PI / 2.0 - 0.01), -(PI / 2.0 - 0.01));

  double r = 1 - (0.0167 * std::cos(0.0172 * (doy - 3)));
  double z = 0.39785 * std::sin(4.868961 + 0.017203 * doy
          + 0.033446 * std::sin(6.224111 + 0.017202 * doy));

  double decl = std::abs(z) < 0.7 ?
      std::atan(z / (std::sqrt(1.0 - z * z))) :
      PI / 2.0 - std::atan(std::sqrt(1 - z * z) / z);

  double z2 = -std::tan(decl) * std::tan(LatRad);

  double h = 0;
  if (z2 < 1.0) {
    if (z2 <= -1.0) {
      h = PI;
    } else {
      double TA = std::abs(z2);
      double AC;
      if (TA < 0.7) {
        AC = 1.570796 - std::atan(TA / std::sqrt(1.0 - TA * TA));
      } else {
        AC = std::atan(std::sqrt(1 - TA * TA) / TA);
      }

      if (z2 < 0) {
        h = 3.141593 - AC;
      } else {
        h = AC;
      }
    }
  }

  double DayLength = 2.0 * (h * 24.0) / (2.0 * PI);
  DayLength = DayLength * 60; // [min]
  return DayLength;
}


double HighTLim(double tleaf) {
  double SHR_CONST_TKFRZ = 273.15;
  double SHR_CONST_RGAS = 8314.467591;

  return 1.0 / (1.0 + std::exp((-2.2 * std::pow(10.0, 5.0) + 710.0 * (tleaf + SHR_CONST_TKFRZ))
          / (SHR_CONST_RGAS * 0.001 * (tleaf + SHR_CONST_TKFRZ))));
}


// This function calculate the net photosynthetic rate based on Farquhar
// model, with updated leaf temperature based on energy balances
void Photosynthesis(double PARi, double LUE, double LER, double pressure, double windv,
                    double tair, double relh, double CO2a, double mp, double Vcmax25,
                    double* A, double* tleaf, double* Resp)
{
  if (tair <= 0.) {
    double ARAD = PARi / 2.3*(1.0 - std::exp(-LER));
    double q10actf = 2.4; // Q10 coefficients

    *tleaf = tair +  ARAD / 38.4;
    *A = 0.;

    double q10act = q10actf * std::exp(-0.009 * (*tleaf - 15.0));
    double Vcmax = Vcmax25 * HighTLim(*tleaf) * std::pow(q10act, 0.1 * (*tleaf - 25.0));
    *Resp = Vcmax * 0.0089; // maintenance respiration
	if (*tleaf < -1.0) (*Resp) *= 0.1;

  } else {
    double ARAD = PARi * (1.0 - std::exp(-LER)) / 2.3;
    double q10actf = 2.4; // Q10 coefficients
    double APAR = PARi * (1.0 - std::exp(-LER)) * 0.95; //assumes only 5% reflectance
    double JmeanL = APAR*LUE*4.0; //4.0 is a factor converting electron to CO2

    double rsmax0 = 2.e4;
    double o2a = 209460.0;
    double co2c = CO2a * pressure * 1.e-6;
    double o2c = o2a * pressure * 1.e-6;
    double phi = 0.98; // coefficient for interpolation, which is normally close to 1.0

    double dleaf = 0.04;
    double kc25 = 30.0;
    double ko25 = 30000.0;
    double akc = 2.1;
    double ako = 1.2;
    double bp = 2000.0;
    double R = 8.314;

    double rb0 = 100.0 * std::sqrt(dleaf / windv);

    double tleafold = tair;

    double tairk = tair + 273.15;
    double es, esdT, qs, qsdT;
    QSat qsat;
    qsat(tairk, pressure, &es, &esdT, &qs, &qsdT);
    double ea = es * relh;

    // output
    double myA;
    double tleafnew = tair;
    double Vcmax;

    // loop to converge tleafnew?
    bool done = false;
    int itr = 0;
    while (!done) {
      itr++;

      tleafold = tleafnew;

      // convert temps to K
      double tleafk = tleafnew + 273.15;

      double q10act = q10actf * std::exp(-0.009 * (tleafnew - 15.0));
      double cf = pressure * 1.e6 / (R * tleafk);
      double rb = rb0 / cf;

      double k_o = ko25 * std::pow(ako, 0.1 * (tleafnew - 25.0));
      double k_c = kc25 * std::pow(akc, 0.1 * (tleafnew - 25.0));
      double c_p = 0.5 * k_c / k_o * o2c * 0.21;
      double awc = k_c * (1.0 + o2c / k_o);

      Vcmax = Vcmax25 * HighTLim(tleafnew) * std::pow(q10act, (0.1 * (tleafnew - 25.0)));
      double We = 0.5 * Vcmax;

      double ei;
      qsat(tleafk, pressure, &ei, &esdT, &qs, &qsdT);
      double cea = std::max(0.3 * ei, std::min(ea, ei));

      // converge ci?
      double ci = 0.7 * co2c;
      double rs;
      for (int i=0; i!=5; ++i) {
        double Kj = std::max(ci - c_p, 0.0) / (4.0 * ci + 8.0 * c_p);
        double Kc = std::max(ci - c_p, 0.0) / (ci + awc);

        double Wc = Kc * Vcmax;
        double Wj = Kj * JmeanL;

        // interpolation for smooth change
        double Wcj = Wc + Wj - std::sqrt((Wc + Wj)*(Wc + Wj) - 4.0*phi*Wc*Wj);
        Wcj = 0.5*Wcj / phi;

        myA = (Wcj + We - std::sqrt((Wcj + We)*(Wcj + We) - 4.0*phi*Wcj*We)) / (2.0*phi);

        // calculate leaf internal [CO2]
        double c_s = 1.0 * std::max(co2c - 1.37 * rb * pressure * myA, 1.e-6);
        double atmp = mp * myA * pressure * cea / (c_s * ei) + bp;
        double btmp = (mp * myA * pressure / c_s + bp) * rb - 1.0;
        double ctmp = -rb;

        double q;
        if (btmp >= 0) {
          q = -0.5 * (btmp + std::sqrt(btmp * btmp - 4.0 * atmp * ctmp));
        } else {
          q = -0.5 * (btmp - std::sqrt(btmp * btmp - 4.0 * atmp * ctmp));
        }
        double r1 = q / atmp;
        double r2 = ctmp / q;
        rs = std::max(r1, r2);
        ci = std::max(c_s - myA * pressure * 1.65 * rs, 0.0);
      }

      double lamda = (2501000 - 2400 * tleafnew) * 14.0 / 1000 * 1.0 / 1000000; // J/kg to J/mol to J/umol;
      tleafnew = tair + 1 / 38.4 * (ARAD - (lamda / rs * (es - ea) / pressure));
      //4.6 is the unit conversion factor from umol photon/m2/s->w/s
      // but there is no 4.6 in this code.........

      // check convergence criteria
      done = itr > 10 || std::abs(tleafnew - tleafold) < 0.1;
    }

    *tleaf = tleafnew;
	*Resp = Vcmax * 0.0089; //maintenance respiration
	if (tleafnew < -1.0) (*Resp) *= 0.1;
    *A = myA;
  }
  return;
}


} // namespace
} // namespace
} // namespace

