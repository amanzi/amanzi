/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#include "EOSFactory.hh"
#include "EOS_VaporInGas.hh"

namespace Amanzi {
namespace AmanziEOS {

Utils::RegisteredFactory<EOS, EOS_VaporInGas> EOS_VaporInGas::factory_("vapor in gas");

}  // namespace AmanziEOS
}  // namespace Amanzi
