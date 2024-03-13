/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Equation of state for viscosity = f(T, P)
*/

#ifndef AMANZI_EOS_VISCOSITY_HH_
#define AMANZI_EOS_VISCOSITY_HH_

namespace Amanzi {
namespace AmanziEOS {

class EOS_Viscosity {
 public:
  EOS_Viscosity(Teuchos::ParameterList& eos_plist) : eos_plist_(eos_plist), ierr_(0){};
  virtual ~EOS_Viscosity(){};

  virtual double Viscosity(double T, double p) = 0;
  virtual double DViscosityDT(double T, double p) = 0;
  virtual double DViscosityDp(double T, double p) = 0;

  // error messages
  int error_code() { return ierr_; }
  std::string error_msg() { return error_msg_; }

 protected:
  Teuchos::ParameterList eos_plist_;

  int ierr_;
  std::string error_msg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
