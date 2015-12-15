/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (berndt@lanl.gov) 
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#include "WRM_vanGenuchten.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM,WRM_vanGenuchten> WRM_vanGenuchten::factory_("van Genuchten");

}  // namespace Flow
}  // namespace Amanzi
