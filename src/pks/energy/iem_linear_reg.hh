/*
  This is the energy component of the ATS and Amanzi codes. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

#include "iem_linear.hh"

namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<IEM,IEMLinear> IEMLinear::factory_("linear");

}  // namespace Energy
}  // namespace Amanzi
