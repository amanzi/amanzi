/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Constant viscosity EOS, defaults to reasonable values for water.
*/

#include "ViscosityConstant.hh"

namespace Amanzi {
namespace AmanziEOS {

ViscosityConstant::ViscosityConstant(Teuchos::ParameterList& eos_plist) : EOS_Viscosity(eos_plist)
{
  InitializeFromPlist_();
};


void
ViscosityConstant::InitializeFromPlist_()
{
  visc_ = eos_plist_.get<double>("viscosity", 8.9e-4);
};

} // namespace AmanziEOS
} // namespace Amanzi
