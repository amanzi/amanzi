/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (berndt@lanl.gov) 
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "WRM_BrooksCorey.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM,WRM_BrooksCorey> WRM_BrooksCorey::factory_("Brooks Corey");

}  // namespace Flow
}  // namespace Amanzi
