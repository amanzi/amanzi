/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for liquid water. See the permafrost physical properties notes for
  references and documentation of this EOS at:
*/

#include "errors.hh"
#include "ViscosityWaterFEHM.hh"

namespace Amanzi {
namespace AmanziEOS {

ViscosityWaterFEHM::ViscosityWaterFEHM(Teuchos::ParameterList& eos_plist)
  : ViscosityBase(eos_plist),
    y0_(0.17409149e-02),  z0_(1.0),
    y1_(0.18894882e-04),  z1_(0.10523153e-01),
    y2_(-0.66439332e-07), z2_(-0.22658391e-05),
    y3_(-0.23122388e-09), z3_(-0.31796607e-06),
    y4_(-0.31534914e-05), z4_(0.29869141e-01),
    y5_(0.11120716e-07),  z5_(0.21844248e-03),
    y6_(-0.48576020e-10), z6_(0.87658855e-06),
    y7_(0.28006861e-07),  z7_(0.41690362e-03),
    y8_(0.23225035e-09),  z8_(-0.25147022e-05),
    y9_(0.47180171e-10),  z9_(0.22144660e-05),
    T0_(273.15) {};


double ViscosityWaterFEHM::Viscosity(double T, double p) {
  double p1(p / 1e+6), t1(T - T0_);
  double p2(p1 * p1), t2(t1 * t1);
  double p3(p2 * p1), t3(t2 * t1);
  double y = y0_ + y1_ * p1 + y2_ * p2 + y3_ * p3 + y4_ * t1 + y5_ * t2 
                 + y6_ * t3 + y7_ * p1 * t1 + y8_ * p2 * t1 + y9_ * p1 * t2;
  double z = z0_ + z1_ * p1 + z2_ * p2 + z3_ * p3 + z4_ * t1 + z5_ * t2 
                 + z6_ * t3 + z7_ * p1 * t1 + z8_ * p2 * t1 + z9_ * p1 * t2;
  return y / z;
}


double ViscosityWaterFEHM::DViscosityDT(double T, double p) {
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
};


double ViscosityWaterFEHM::DViscosityDp(double T, double p) {
  Errors::Message message("EOS viscosity of water: derivative not implemented");
  Exceptions::amanzi_throw(message);
  return -1.0;
};

}  // namespace AmanziEOS
}  // namespace Amanzi
