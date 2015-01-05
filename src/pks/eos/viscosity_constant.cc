/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Constant viscosity EOS, defaults to reasonable values for water.
*/

#include "viscosity_constant.hh"

namespace Amanzi {
namespace Relations {

ViscosityConstant::ViscosityConstant(Teuchos::ParameterList& visc_plist) :
    visc_plist_(visc_plist) {
  InitializeFromPlist_();
};


void ViscosityConstant::InitializeFromPlist_() {
  // defaults to water
  visc_ = visc_plist_.get<double>("Viscosity [kg/m-s]", 8.9e-4);
};

}  // namespace Relations
}  // namespace Amanzi
