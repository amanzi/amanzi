/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Constant density/viscosity EOS, defaults to reasonable values for water.
*/

#include "eos_constant.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<EOS,EOSConstant> EOSConstant::factory_("constant");

}  // namespace Relations
}  // namespace Amanzi
