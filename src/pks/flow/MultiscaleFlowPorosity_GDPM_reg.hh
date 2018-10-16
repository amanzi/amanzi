/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#include "MultiscaleFlowPorosity_GDPM.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<MultiscaleFlowPorosity, MultiscaleFlowPorosity_GDPM>
    MultiscaleFlowPorosity_GDPM::factory_("generalized dual porosity");

}  // namespace Flow
}  // namespace Amanzi
