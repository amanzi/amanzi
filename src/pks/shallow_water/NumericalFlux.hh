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
  double minmod(double a, double b);

  virtual std::vector<double> Compute(
          const std::vector<double>& UL, const std::vector<double>& UR) = 0;

  double MaxSpeed() const { return lambda_max_; }
  double MinSpeed() const { return lambda_min_; }

 protected:
  double g_;
  double lambda_max_, lambda_min_;
  int shallow_water_model_;
  double pipe_diameter_;
  double pipe_cross_section_;
  double celerity_;

};


inline
std::vector<double> NumericalFlux::PhysicalFlux(const std::vector<double>& U)
{
  std::vector<double> F(3);

  double h, h2, u, v, WettedAngle;
  double eps2 = 1e-12;

  // transform from conservative (h, hu, hv) to primary (h, u, v) variables
  h = U[0];
  h2 = h * h;
  u = 2.0 * h * U[1] / (h2 + std::fmax(h2, eps2));
  v = 2.0 * h * U[2] / (h2 + std::fmax(h2, eps2));
  // get reconstructed wetted angle
  WettedAngle = U[3];

  double HydrostaticPressureForce = 0.0;

  if (shallow_water_model_){

     HydrostaticPressureForce = 0.5 * g_ * h2;

  }

  else {

     if ((0.0 < h && h < pipe_cross_section_)){ //flow is ventilated (free-surface)

        HydrostaticPressureForce = 3.0 * sin(WettedAngle * 0.5) - pow(sin(WettedAngle * 0.5),3) - 3.0 * (WettedAngle * 0.5) * cos(WettedAngle * 0.5);
        HydrostaticPressureForce = HydrostaticPressureForce * g_ * pow(pipe_diameter_,3) / 24.0;

     }

     else if (h >= pipe_cross_section_) { //flow is pressurized

         double PressurizedHead = (celerity_ * celerity_ * (h - pipe_cross_section_)) / (g_ * pipe_cross_section_);
         HydrostaticPressureForce = g_ * h * (PressurizedHead + sqrt(h/3.14159265359)); 

      }

      else if (h < 0.0) {
          std::cout  << "negative wetted area in NumericalFlux.h " << std::endl;
          abort();
      }

  }

  F[0] = h * u;
  F[1] = h * u * u + HydrostaticPressureForce;
  F[2] = h * u * v;

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
  
