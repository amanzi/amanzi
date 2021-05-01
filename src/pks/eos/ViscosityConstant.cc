/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Constant viscosity EOS, defaults to reasonable values for water.
*/

#include "ViscosityConstant.hh"

namespace Amanzi {
namespace AmanziEOS {

ViscosityConstant::ViscosityConstant(Teuchos::ParameterList& eos_plist)
  : ViscosityBase(eos_plist) {
  InitializeFromPlist_();
};


void ViscosityConstant::InitializeFromPlist_() {
  visc_ = eos_plist_.get<double>("viscosity", 8.9e-4);
};

}  // namespace AmanziEOS
}  // namespace Amanzi
