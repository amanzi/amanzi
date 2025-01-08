/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  EOS

  EOS is purely virtual base class for an equation of state.
*/

#ifndef AMANZI_EOS_DENSITY_HH_
#define AMANZI_EOS_DENSITY_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace AmanziEOS {

class EOS_Density {
 public:
  EOS_Density(Teuchos::ParameterList& eos_plist) : ierr_(0)
  {
    M_ = eos_plist.get<double>("molar mass");

    if (eos_plist.isParameter("molar density"))
      rho_ = eos_plist.get<double>("molar density") * M_;
    else
      rho_ = eos_plist.get<double>("density");
  }
  virtual ~EOS_Density(){};

  // Virtual methods that form the EOS
  virtual double Density(double T, double p) = 0;
  virtual double DDensityDT(double T, double p) = 0;
  virtual double DDensityDp(double T, double p) = 0;

  virtual double MolarDensity(double T, double p) = 0;
  virtual double DMolarDensityDT(double T, double p) = 0;
  virtual double DMolarDensityDp(double T, double p) = 0;

  double MolarMass() { return M_; }

  // error messages
  int error_code() { return ierr_; }
  std::string error_msg() { return error_msg_; }

 protected:
  double M_, rho_;

  int ierr_;
  std::string error_msg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
