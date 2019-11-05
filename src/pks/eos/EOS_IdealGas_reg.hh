/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for the ideal gas. It does not implement viscosity at this point!
*/

#include "EOS_IdealGas.hh"

namespace Amanzi {
namespace AmanziEOS {

// registry of method
Utils::RegisteredFactory<EOS, EOS_IdealGas> EOS_IdealGas::factory_("ideal gas");

}  // namespace AmanziEOS
}  // namespace Amanzi
