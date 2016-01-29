/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for TC implementations.
*/

#include "TCMFactory_TwoPhase.hh"

// explicity instantitate the static data of Factory<TCM_TwoPhase>
namespace Amanzi {
namespace Utils {

template<> 
Factory<Energy::TCM_TwoPhase>::map_type* Factory<Energy::TCM_TwoPhase>::map_;

}  // namespace Utils
}  // namespace Amanzi

