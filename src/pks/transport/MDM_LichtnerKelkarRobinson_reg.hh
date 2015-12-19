/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "MDM_LichtnerKelkarRobinson.hh"

namespace Amanzi {
namespace Transport {

Utils::RegisteredFactory<MDM,MDM_LichtnerKelkarRobinson> 
    MDM_LichtnerKelkarRobinson::factory_("Lichtner-Kelkar-Robinson");

}  // namespace Transport
}  // namespace Amanzi
