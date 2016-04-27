/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Constant density/viscosity EOS. It defaults to reasonable values 
  for water.
*/

#include "EOS_Constant.hh"

namespace Amanzi {
namespace EOS {

// registry of method
Utils::RegisteredFactory<EOS, EOS_Constant> EOS_Constant::factory_("constant");

}  // namespace EOS
}  // namespace Amanzi
