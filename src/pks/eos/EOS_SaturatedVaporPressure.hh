/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Base class for saturated vapor pressure.
*/

#ifndef AMANZI_EOS_SATURATED_VAPOR_PRESSURE_HH_
#define AMANZI_EOS_SATURATED_VAPOR_PRESSURE_HH_

namespace Amanzi {
namespace AmanziEOS {

class EOS_SaturatedVaporPressure {
 public:
  virtual ~EOS_SaturatedVaporPressure() {};

  virtual double Pressure(double T) = 0;
  virtual double DPressureDT(double T) = 0;

  // error messages
  int error_code() { return ierr_; }
  std::string error_msg() { return error_msg_; }

 protected:
  int ierr_;
  std::string error_msg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
