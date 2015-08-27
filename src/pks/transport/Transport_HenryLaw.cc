/*
  This is the transport component of Amanzi. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//TPLs
#include <boost/algorithm/string.hpp>

#include "Transport_PK.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
* Re-partition components between liquid and gas phases.
******************************************************************* */
void Transport_PK::MakeAirWaterPartitioning_()
{
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell", false);
  int num_components = tcc_next.NumVectors();

  for (int ig = num_aqueous; ig < num_components; ig++) {
    std::string name_l = boost::replace_all_copy(component_names_[ig], "(g)", "(l)");
    int il = FindComponentNumber(name_l);
    if (il < 0 || il >= num_aqueous) {
      Errors::Message msg;
      msg << "Gas component \"" << component_names_[ig] 
          << "\" has no matching liquid component \"" << name_l << "\"\n";
      Exceptions::amanzi_throw(msg);
    }
  }
}

}  // namespace Transport
}  // namespace Amanzi


