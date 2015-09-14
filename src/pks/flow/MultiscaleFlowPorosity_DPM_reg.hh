/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#include "MultiscaleFlowPorosity_DPM.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<MultiscaleFlowPorosity, MultiscaleFlowPorosity_DPM>
    MultiscaleFlowPorosity_DPM::factory_("dual porosity");

}  // namespace Flow
}  // namespace Amanzi
