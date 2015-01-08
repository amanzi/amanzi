/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  EOS for an ideal gas (does not implement viscosity at this point!)
*/

#include "eos_factory.hh"
#include "eos_vapor_in_gas.hh"

namespace Amanzi {
namespace Relations {

Utils::RegisteredFactory<EOS,EOSVaporInGas> EOSVaporInGas::factory_("vapor in gas");

}  // namespace Relations
}  // namespace Amanzi
