/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Registry of a viscosity method.
*/

#include "Viscosity_Constant.hh"

namespace Amanzi {
namespace AmanziEOS {

Utils::RegisteredFactory<Viscosity_Base, Viscosity_Constant> Viscosity_Constant::factory_("constant");

}  // namespace AmanziEOS
}  // namespace Amanzi
