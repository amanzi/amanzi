/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

#include "IEM_Linear.hh"

namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<IEM,IEM_Linear> IEM_Linear::factory_("linear");

}  // namespace Energy
}  // namespace Amanzi
