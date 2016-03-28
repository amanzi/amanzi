/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#include "MultiscaleTransportPorosity_DPM.hh"

namespace Amanzi {
namespace Transport {

Utils::RegisteredFactory<MultiscaleTransportPorosity, MultiscaleTransportPorosity_DPM>
    MultiscaleTransportPorosity_DPM::factory_("dual porosity");

}  // namespace Transport
}  // namespace Amanzi
