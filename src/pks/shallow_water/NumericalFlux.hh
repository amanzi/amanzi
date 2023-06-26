/*
  Shallow Water PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef SHALLOW_WATER_NUMERICAL_FLUX_HH_
#define SHALLOW_WATER_NUMERICAL_FLUX_HH_

#include <vector>

namespace Amanzi {
namespace ShallowWater {

class NumericalFlux {
 public:
  virtual ~NumericalFlux() {};

  std::vector<double> PhysicalFlux(const std::vector<double>& U);
  double HydrostaticPressureForce(const double & PrimaryVariable, const double & WettedAngle);
  double minmod(double a, double b);

  virtual std::vector<double> Compute(
          const std::vector<double>& UL, const std::vector<double>& UR) = 0;

  double MaxSpeed() const { return lambda_max_; }
  double MinSpeed() const { return lambda_min_; }

 protected:
  double g_;
  double lambda_max_, lambda_min_;
  double pipe_diameter_;
  double pipe_cross_section_;
  double celerity_;

};

inline 
double NumericalFlux::HydrostaticPressureForce(const double & PrimaryVariable, const double & WettedAngle)
{

  double HydroPressForce = 0.0;

  if (WettedAngle < 0.0){ // shallow water 

     HydroPressForce = 0.5 * g_ * PrimaryVariable * PrimaryVariable;

  }

  else { // pipe flow

     if ((0.0 < PrimaryVariable && PrimaryVariable < pipe_cross_section_)){ //flow is ventilated (free-surface)

        HydroPressForce = 3.0 * sin(WettedAngle * 0.5) - pow(sin(WettedAngle * 0.5),3) - 3.0 * (WettedAngle * 0.5) * cos(WettedAngle * 0.5);
        HydroPressForce = HydroPressForce * g_ * pow(pipe_diameter_,3) / 24.0;

     }

     else if (PrimaryVariable >= pipe_cross_section_) { //flow is pressurized

         double PressureHead = (celerity_ * celerity_ * (PrimaryVariable - pipe_cross_section_)) / (g_ * pipe_cross_section_);
         HydroPressForce = g_ * PrimaryVariable * (PressureHead + sqrt(PrimaryVariable / 3.14159265359));

      }

      else if (PrimaryVariable < 0.0) {
          std::cout  << "negative wetted area in NumericalFlux.h (HydrostaticPressureForce) " << std::endl;
          abort();
      }

  }

  return HydroPressForce;

}

inline
std::vector<double> NumericalFlux::PhysicalFlux(const std::vector<double>& U)
{
  std::vector<double> F(3);

  double h2, u;
  double eps2 = 1e-12;

  // U[0] = primary variable (ponded depth or wetted area)
  // U[1] = U[0] * u
  // U[2] = U[0] * v
  // U[3] = wetted angle
  h2 = U[0] * U[0];
  u = 2.0 * U[0] * U[1] / (h2 + std::fmax(h2, eps2));

  // hydrostatic pressure force
  double HPF = HydrostaticPressureForce(U[0], U[3]);

  F[0] = U[1];
  F[1] = U[1] * u + HPF;
  F[2] = U[2] * u;

  return F;
}


//--------------------------------------------------------------
// minmod function
//--------------------------------------------------------------
inline
double NumericalFlux::minmod(double a, double b)
{
  double m;

  if (a*b > 0) {
    if (std::fabs(a) < std::fabs(b)) {
      m = a;
    } else {
      m = b;
    }
  } else {
    m = 0.0;
  }

  return m;
}


}  // namespace ShallowWater
}  // namespace Amanzi
  
#endif
  
