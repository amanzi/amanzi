/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Equation of state for diffusion coefficient d(T, P)
*/

#ifndef AMANZI_EOS_DIFFUSION_HH_
#define AMANZI_EOS_DIFFUSION_HH_

namespace Amanzi {
namespace AmanziEOS {

class EOS_Diffusion {
 public:
  EOS_Diffusion(Teuchos::ParameterList& plist) : plist_(plist), ierr_(0){};
  virtual ~EOS_Diffusion(){};

  virtual double Diffusion(double T, double p) = 0;
  virtual double DDiffusionDT(double T, double p) = 0;
  virtual double DDiffusionDp(double T, double p) = 0;

  // error messages
  int error_code() { return ierr_; }
  std::string error_msg() { return error_msg_; }

 protected:
  Teuchos::ParameterList plist_;

  int ierr_;
  std::string error_msg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
