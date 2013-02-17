/*
This is the transport component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors:  Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "Transport_BC_Factory.hh"
#include "Transport_constants.hh"

namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
* Process Dirichet BC (concentration), step 1.
****************************************************************** */
void TransportBCFactory::CreateConcentration(
    std::vector<BoundaryFunction*>& bcs, std::vector<int> bcs_tcc_index) const
{
}

}  // namespace AmanziTransport
}  // namespace Amanzi
