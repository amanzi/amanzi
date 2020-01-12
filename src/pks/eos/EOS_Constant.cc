/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Constant density/viscosity EOS. It defaults to reasonable values 
  for water at 25C.
*/

#include "EOS_Constant.hh"

namespace Amanzi {
namespace AmanziEOS {

EOS_Constant::EOS_Constant(Teuchos::ParameterList& eos_plist) :
    eos_plist_(eos_plist) {
  InitializeFromPlist_();
};


void EOS_Constant::InitializeFromPlist_() {
  // defaults to water
  M_ = eos_plist_.get<double>("molar mass", 18.0153e-03);

  if (eos_plist_.isParameter("molar density")) {
    rho_ = eos_plist_.get<double>("molar density") * M_;
  } else {
    rho_ = eos_plist_.get<double>("density", 997.07);
  }
};

}  // namespace AmanziEOS
}  // namespace Amanzi
