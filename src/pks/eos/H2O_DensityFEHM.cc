/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for liquid water for T between 0.001 and 360 C from FEHM manual
*/

#include "H2O_DensityFEHM.hh"

namespace Amanzi {
namespace AmanziEOS {

H2O_DensityFEHM::H2O_DensityFEHM(Teuchos::ParameterList& eos_plist)
  : EOS_Density(eos_plist),
    y0_(1.0),            z0_(1.0009476e-03),
    y1_(1.7472599e-02),  z1_(1.6812589e-05),
    y2_(-2.0443098e-05), z2_(-2.4582622e-08),
    y3_(1.7442012e-07),  z3_(-1.7014984e-10),
    y4_(4.9564109e-03),  z4_(4.8841156e-06),
    y5_(-4.0757664e-05), z5_(-3.2967985e-08),
    y6_(5.0676664e-08),  z6_(2.861938e-11),
    y7_(5.0330978e-05),  z7_(5.3249055e-08),
    y8_(3.3914814e-07),  z8_(3.0456698e-10),
    y9_(-1.8383009e-07), z9_(-1.2221899e-10),
    T0_(273.15) {};


double H2O_DensityFEHM::Density(double T, double p) {
  double p1(p / 1e+6), t1(T - T0_);
  double p2(p1 * p1), t2(t1 * t1);
  double p3(p2 * p1), t3(t2 * t1);
  double y = y0_ + y1_ * p1 + y2_ * p2 + y3_ * p3 + y4_ * t1 + y5_ * t2 
                 + y6_ * t3 + y7_ * p1 * t1 + y8_ * p2 * t1 + y9_ * p1 * t2;
  double z = z0_ + z1_ * p1 + z2_ * p2 + z3_ * p3 + z4_ * t1 + z5_ * t2 
                 + z6_ * t3 + z7_ * p1 * t1 + z8_ * p2 * t1 + z9_ * p1 * t2;
  return y / z;
}


double H2O_DensityFEHM::DDensityDT(double T, double p) {
  double p1(p / 1e+6), t1(T - T0_);
  double p2(p1 * p1), t2(t1 * t1);
  double p3(p2 * p1), t3(t2 * t1);

  double dydt = y4_ + 2 * y5_ * t1 + 3 * y6_ * t2 + y7_ * p1 + y8_ * p2 + 2 * y9_ * p1 * t1;
  double dzdt = z4_ + 2 * z5_ * t1 + 3 * z6_ * t2 + z7_ * p1 + z8_ * p2 + 2 * z9_ * p1 * t1;

  double y = y0_ + y1_ * p1 + y2_ * p2 + y3_ * p3 + y4_ * t1 + y5_ * t2 
                 + y6_ * t3 + y7_ * p1 * t1 + y8_ * p2 * t1 + y9_ * p1 * t2;
  double z = z0_ + z1_ * p1 + z2_ * p2 + z3_ * p3 + z4_ * t1 + z5_ * t2 
                 + z6_ * t3 + z7_ * p1 * t1 + z8_ * p2 * t1 + z9_ * p1 * t2;
  return (dydt * z - dzdt * y) / (z * z);
}


double H2O_DensityFEHM::DDensityDp(double T, double p) {
  double p1(p / 1e+6), t1(T - T0_);
  double p2(p1 * p1), t2(t1 * t1);
  double p3(p2 * p1), t3(t2 * t1);

  double dydp = y1_ + 2 * y2_ * p1 + 3 * y3_ * p2 + y7_ * t1 + 2 * y8_ * p1 * t1 + y9_ * t2;
  double dzdp = z1_ + 2 * z2_ * p1 + 3 * z3_ * p2 + z7_ * t1 + 2 * z8_ * p1 * t1 + z9_ * t2;

  double y = y0_ + y1_ * p1 + y2_ * p2 + y3_ * p3 + y4_ * t1 + y5_ * t2 
                 + y6_ * t3 + y7_ * p1 * t1 + y8_ * p2 * t1 + y9_ * p1 * t2;
  double z = z0_ + z1_ * p1 + z2_ * p2 + z3_ * p3 + z4_ * t1 + z5_ * t2 
                 + z6_ * t3 + z7_ * p1 * t1 + z8_ * p2 * t1 + z9_ * p1 * t2;
  return (dydp * z - dzdp * y) / (z * z) / 1e+6;
}

}  // namespace AmanziEOS
}  // namespace Amanzi
