/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "errors.hh"
#include "tabular-function.hh"

#include "Flow_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ****************************************************************
* Process string for the discretization method.
**************************************************************** */
void Flow_PK::ProcessStringMFD3D(const std::string name, int* method)
{
  if (name == "monotone") {
    *method = FLOW_MFD3D_HEXAHEDRA_MONOTONE;
  } else if (name == "support operator") {
    *method = FLOW_MFD3D_SUPPORT_OPERATOR;
  } else if (name == "two point flux approximation") {
    *method = FLOW_MFD3D_TWO_POINT_FLUX;
  } else if (name == "optimized") {
    *method = FLOW_MFD3D_OPTIMIZED;
  } else {
    *method = FLOW_MFD3D_POLYHEDRA;
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

